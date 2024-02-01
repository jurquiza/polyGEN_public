from flask import Flask, redirect, url_for, render_template, request, Response, send_from_directory
from io import BytesIO
from zipfile import ZipFile
from datetime import date,time
import os, tempfile

from engine import *
    
session = {}		# dictionary to pass variables around routes


app = Flask(__name__)
app.secret_key = 'ee839687a282e6493054c86e00e86925dc04de931acf4aed'

@app.route("/")
def home():
    return render_template("main.html")


@app.route("/learn")
def learn():
    return render_template("learn_more.html")


@app.route("/ptg", methods=["POST","GET"])
def sequence():

    ## Setting up defaults and variables
    session['clr'] = {'sequence_spacers': '#FFFFFF', 'link': '#FFFFFF', 'poltype_input': '#FFFFFF', 'oligo_index': '#FFFFFF', 'PTG_name': '#FFFFFF'}
    enzms={'bsai': ['gaggtctcg', 'cgagacctc'], 'bsmbi': ['tgcgtctca', 'tgagacgca'], 'btgzi': ['ctgcgatggagtatgtta', 'taacatactccatcgcag'], 'bbsi': ['ttgaagactt', 'aagtcttcaa']} #templates found in pUU080 (bsai), pUPD2 (bsmbi), Ortega-Escalante et al. 2018 (btgzi), Kun (bbsi)
    
    if request.method == "POST":
        
        if request.form['submitPTG'] == 'submit':
            
            ## Pulling all input from the page
            session['msg'] = None
            session['poltype'] = request.form["poltype_input"]
            session['enzm'] = request.form['enzm_input']
            session['PTG_name'] = request.form["PTG_name"].replace(" ", "_")
            session['oligo_prefix'] = request.form["oligo_prefix"]
            session['oligo_index'] = request.form["oligo_index"]
            session['PTG_transfer'] = request.form["sequence_spacers"]
            session['bb_linkers'] = request.form["bb_linkers"]
            session['ad_linkers'] = request.form["ad_linkers"]
            session['tm_range'] = [int(request.form["min_temp"][:2]), int(request.form["max_temp"][:2])]
            session['staticBorderPrimers'] = request.form.get('staticBorderPrimers') if request.form.get('staticBorderPrimers') else False
            session['noBorderPrimers'] = request.form.get('noBorderPrimers') if request.form.get('noBorderPrimers') else False
            
            if not session['PTG_name']:
                session['PTG_name'] = request.form['poltype_input'].upper()
            if not session['oligo_prefix']:
                session['oligo_prefix'] = 'o'
            if not session['oligo_index']:
                session['oligo_index'] = '0'
            
            
            ## Catching input errors
            if '/' in session['PTG_name']:
                raise InvalidUsage("Polycistron name may not contain a '/'", status_code=400, payload={'pge': 'sequence.html', 'box': 'PTG_name'})
            if session['PTG_transfer'] == '':
                raise InvalidUsage("You must specify a polycistron description", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
            if session['oligo_index'] != '' and re.search(r'^[0-9]*$', session['oligo_index']) is None:
                raise InvalidUsage("Starting index must be a number", status_code=400, payload={'pge': 'sequence.html', 'box': 'oligo_index'})
            
            bbLinkersSplit = session['bb_linkers'].split(';') if session['bb_linkers'] else ['tgcc', 'gttt']
            adLinkersSplit = session['ad_linkers'].split(';') if session['ad_linkers'] else []
            for lnk in bbLinkersSplit+adLinkersSplit:
                if len(lnk) != 4 or re.search(r'^[ACGTacgt]*$', lnk) is None:
                    raise InvalidUsage("Invalid linker input", status_code=400, payload={'pge': 'sequence.html', 'box': 'link'})
            
            
            ## Preprocessing input and catching input errors. The error catching might be duplicated in the individual functions
            ## but catching the errors earlier saves computing time.
            PTG_input = session['PTG_transfer'].replace(" ", "").replace("\r\n", "").split('|')
            PTG_structure = []
            for element in PTG_input:
                e = element.split(';')
                if len(e) != 2:
                    raise InvalidUsage("Invalid input syntax", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
                elif e[0] not in ['gRNA', 'pegRNA', 'smRNA', 'crRNA']:
                    raise InvalidUsage("Invalid RNA type", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
                elif re.search(r'^[ACGTacgt]*$', e[1]) is None:
                    raise InvalidUsage("Invalid sequence input", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
                else:
                    PTG_structure.append(e)
            
            
            ## Running inputs through engine functions   
            PTG = PTGbldr(session['PTG_name'], PTG_structure, session['poltype'])
            
            ggArgs = {'poltype': session['poltype'],
                      'enzm': session['enzm'],
                      'tm_range': session['tm_range'],
                      'bb_linkers': bbLinkersSplit,
                      'ad_linkers': adLinkersSplit}
            polycistron = scarless_gg(PTG, **ggArgs)
            
            primerArgs = {'oligo_prefix': session['oligo_prefix'], 
                          'oligo_index': session['oligo_index'], 
                          'staticBorderPrimers': session['staticBorderPrimers'],
                          'noBorderPrimers': session['noBorderPrimers'], 
                          'poltype': session['poltype'], 
                          'enzm': session['enzm'],
                          'bb_linkers': ggArgs['bb_linkers'],
                          'ad_linkers': ggArgs['ad_linkers']}
            session['plcstrn'] = annotatePrimers(polycistron, **primerArgs)
            
            if session['plcstrn'].warning:
                session['msg'] = session['plcstrn'].warning
            
            return render_template("primer_list.html", session=session)
            
            
        ## If the reset button is pressed, erase all input
        elif request.form['submitPTG'] == 'reset':
            session['PTG_name'] = session['oligo_prefix'] = session['oligo_index'] = session['bb_linkers'] = session['ad_linkers'] = session['PTG_transfer'] = ''
            return render_template("sequence.html",PTG_transfer=session.get('PTG_transfer', None), session=session)
            
    else:
        
        return render_template("sequence.html",PTG_transfer=session.get('PTG_transfer', None), session=session)


@app.route("/peg", methods=["POST","GET"])
def peg_generation():
    
    ## Setting up defaults and variables
    session['msg'] = None
    session['clr'] = {'sequence': '#FFFFFF', 'edits': '#FFFFFF'}

    if request.method == "POST":

        if request.form['submitPEG'] == 'submit':
        
            ## Pulling all inputs from page
            session['PEG_sequence'] = request.form["sequence"]
            session['PEG_edits'] = request.form["edits"]
            PEG_mode = request.form["mode"]

            ## Preprocessing input
            # Clean up whitespaces and newlines in sequence
            session['PEG_sequence'] = session['PEG_sequence'].replace(' ', '')
            session['PEG_sequence'] = session['PEG_sequence'].replace('\r\n', '')

            # Concatenate edits
            PEG_edits = session['PEG_edits'].split('|')
            pegs_list = []
            for edt in PEG_edits:
                pegs_list.append(edt.split(';'))

            ## Running inputs through engine function
            edits_list,session['msg'] = pegbldr(session['PEG_sequence'], pegs_list, PEG_mode)

            ## Postprocessing output
            session['PTG_transfer'] = ''
            for c,peg in enumerate(edits_list):
                session['PTG_transfer'] += str(peg[1]) + ';' +  str(peg[2])
                if c < len(edits_list)-1:
                    session['PTG_transfer'] += '|'

            return redirect(url_for('sequence'))

        ## If the reset button is pressed, erase all input
        elif request.form['submitPEG'] == 'reset':
            session['PEG_sequence'] = session['PEG_edits'] = ''
            return render_template("peg_generation.html", PEG_transfer=session.get('PEG_sequence', None), session=session)

    else:
        
        return render_template("peg_generation.html", PEG_transfer=session.get('PEG_sequence', None), session=session)
      
        
@app.route("/primer_list", methods=["POST","GET"])
def serve_primers():
    
    ## Starting to write the csv file
    csv = 'List of oligos:\n'
    csv += 'oligo_ID,sequence\n'
    
    ## Appending the oligo list to the csv file
    for oligo in session['plcstrn'].oligos:
        csv += oligo[0] + ',' + oligo[1] + '\n'
        
    ## Appending the fragment table to the csv file
    csv += '\nTable of fragments:\n'
    csv += 'fragment_id,fragment_type,forward_primer,Tm_forw,reverse_primer,Tm_rev\n'
    for c,fragment in enumerate(session['plcstrn'].parts):
        csv += session['PTG_name']+'_f'+str(c)+','+fragment.type+','+session['plcstrn'].oligos[c*2][0]+','+str(np.round(fragment.primer_forward_tm,1))+','+session['plcstrn'].oligos[c*2+1][0]+','+str(np.round(fragment.primer_reverse_tm, 1))+'\n'
    
    ## Appending the input parameters to the csv file
    csv += '\nInput parameters:\n'
    csv += 'Polycistron name:,' + session['PTG_name'] + '\n'
    csv += 'Oligos prefix:,' + session['oligo_prefix'] + '\n'
    csv += 'Starting index:,' + session['oligo_index'] + '\n'
    csv += 'Border linkers:,' + '+'.join(session['bb_linkers'] if session['bb_linkers'] else ['tgcc', 'gttt']) + '\n'
    csv += 'Addidional linkers:,' + '+'.join(session['ad_linkers'] if session['ad_linkers'] else []) + '\n'
    csv += 'Type of Polycistron:,' + session['poltype'] + '\n'
    csv += 'Restriction enzyme:,' + session['enzm'] + '\n'
    csv += 'Melting temperature range:,' + str(session['tm_range'][0]) + '-' + str(session['tm_range'][1]) + '\n'
    csv += 'Invariable border primers:,' + str(session['staticBorderPrimers']) + '\n'
    csv += 'Omit border primers:,' + str(session['noBorderPrimers']) + '\n'
    
    csv += '\nSequence input:\n'
    PTG_input = session['PTG_transfer'].replace(" ", "").replace("\r\n", "").split("|")
    for part in PTG_input:
        csv += part.replace(";", ",") + '\n'
    
    ## Generating the object for the genbank file
    sr = SeqRecord(seq=Seq(session['plcstrn'].sequence, alphabet=IUPAC.ambiguous_dna), name=session['PTG_name'], annotations={'date': date.today().strftime("%d-%b-%Y").upper(), 'topology': 'linear'})
    for ftr in session['plcstrn'].features:
        sr.features.append(ftr)
    
    ## Jsonifying the raw output
    gb_json = {}
    gb_json['polycistron'] = polyToJson(session['plcstrn'])
    gb_json['msg'] = session['msg']
    
    ## Writing everything to a zip file
    in_memory = BytesIO()
    zf = ZipFile(in_memory, mode='w')
    zf.writestr(session['PTG_name']+"_oligos.csv", csv)
    zf.writestr(session['PTG_name']+".gb", sr.format('genbank'))
    zf.writestr(session['PTG_name']+'_raw.json', json.dumps(gb_json))
    with open('protocol.txt') as f:
        zf.writestr('protocol.txt', f.read())
    zf.close()
    in_memory.seek(0)
    outpt = in_memory.read()
    
    ## Returning the zip file
    return Response(
        outpt,
        mimetype="application/zip",
        headers={'Content-Disposition': 'attachment;filename=pg_%s_%s.zip' %(session['PTG_name'], str(date.today()))})


@app.route("/impressum")
def impress():
    return render_template("impressum.html")


@app.route("/success")
def success():
    return "Success"


@app.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    session['msg'] = error.to_dict()['message']
    session['clr'][error.to_dict()['box']] = '#FA5858'
    return render_template(error.to_dict()['pge'], PEG_transfer=session.get('PEG_sequence'), PTG_transfer=session.get('PTG_transfer', None), session=session)


if __name__=='__main__':
    app.run(debug=True, host='0.0.0.0', port=5000) ##host 0.0.0.0 is important for docker container
