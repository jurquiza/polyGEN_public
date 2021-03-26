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
    session['msg'] = None
    session['clr'] = {'sequence_spacers': '#FFFFFF', 'link': '#FFFFFF', 'poltype_input': '#FFFFFF', 'oligo_index': '#FFFFFF', 'PTG_name': '#FFFFFF'}
    session['enzm_site'] = ['gaggtctcg', 'cgagacctc']
    enzms={'bsai': ['gaggtctcg', 'cgagacctc'], 'bsmbi': ['tgcgtctca', 'tgagacgca'], 'btgzi': ['ctgcgatggagtatgtta', 'taacatactccatcgcag'], 'bbsi': ['agaagacag', 'ctgtcttct']} #templates found in pUU080 (bsai), pUPD2 (bsmbi), Ortega-Escalante et al. 2018 (btgzi), pUU256 (bbsi)
    
    if request.method == "POST":
        if request.form['submit_button'] == 'submit':
            args = {}
            args['poltype'] = request.form["poltype_input"]
            args['enzm'] = request.form['enzm_input']
            session['enzm_site'] = enzms[request.form['enzm_input']]
            session['PTG_name'] = request.form["PTG_name"]
            session['oligo_prefix'] = request.form["oligo_prefix"]
            session['oligo_index'] = request.form["oligo_index"]
            session['PTG_transfer'] = request.form["sequence_spacers"]
            session['bb_ovrhng'] = request.form["bb_ovrhng"]
            session['add_ovrhng'] = request.form["add_ovrhng"]
            session['max_ann_len'] = request.form['maxAnnLen']
        
            args['tm_range'] = [int(request.form["min_temp"][:2]), int(request.form["max_temp"][:2])]
            args['bb_overlaps'] = ['tgcc', 'gttt']
            args['additional_overhangs'] = []
            if session['bb_ovrhng']:
                args['bb_overlaps'] = session["bb_ovrhng"].split(';')
            if session['add_ovrhng']:
                args['additional_overhangs'] = session["add_ovrhng"].split(';')
            if session['max_ann_len']:
                args['max_ann_len'] = int(session['max_ann_len'])
            
            if '/' in session['PTG_name']:
                raise InvalidUsage("Polycistron name may not contain a '/'", status_code=400, payload={'pge': 'sequence.html', 'box': 'PTG_name'})
            if session['PTG_transfer'] == '':
                raise InvalidUsage("You must specify a polycistron description", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
            if session['oligo_index'] != '' and re.search(r'^[0-9]*$', session['oligo_index']) is None:
                raise InvalidUsage("Starting index must be a number", status_code=400, payload={'pge': 'sequence.html', 'box': 'oligo_index'})
            
            PTG_input = session['PTG_transfer'].split('|')
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
            
            for lnk in args['bb_overlaps']+args['additional_overhangs']:
                print(lnk)
                if len(lnk) != 4 or re.search(r'^[ACGTacgt]*$', lnk) is None:
                    raise InvalidUsage("Invalid linker input", status_code=400, payload={'pge': 'sequence.html', 'box': 'link'})
            
            PTG = PTGbldr(session['PTG_name'], PTG_structure, args['poltype'])
            
            print([prt.to_json() for prt in PTG])
            
            session['plcstrn'],session['msg'] = scarless_gg(PTG, **args)
        
            if not session['PTG_name']:
                session['PTG_name'] = request.form['poltype_input'].upper()
            if not session['oligo_prefix']:
                session['oligo_prefix'] = 'o'
            if not session['oligo_index']:
                session['oligo_index'] = '0'
            
            if request.form.get('borderPrimers'):
                if args['poltype_run'] == 'ptg':
                    session['plcstrn'].parts[0].primer_forward = session['enzm_site'][0] + reverse_complement(args['bb_overlaps'][0].upper()) + 'aacaaagcaccagtggtctagtggtag'
                    session['plcstrn'].parts[-1].primer_reverse = reverse_complement(session['enzm_site'][1]) + reverse_complement(args['bb_overlaps'][1].upper()) + 'tgcaccagccgggaatcgaac'
                if args['poltype_run'] == 'ca':
                    session['plcstrn'].parts[0].primer_forward = session['enzm_site'][0] + reverse_complement(args['bb_overlaps'][0].upper()) + 'aatttctactgttgtagat'
                    session['plcstrn'].parts[-1].primer_reverse = reverse_complement(session['enzm_site'][1]) + reverse_complement(args['bb_overlaps'][1].upper()) + 'atctacaacagtagaaatt'

            return render_template("primer_list.html", session=session)
            
        elif request.form['submit_button'] == 'reset':
            session['PTG_name'] = session['oligo_prefix'] = session['oligo_index'] = session['bb_ovrhng'] = session['add_ovrhng'] = session['PTG_transfer'] = ''
            return render_template("sequence.html",PTG_transfer=session.get('PTG_transfer', None), session=session)
            
    else:
        
        return render_template("sequence.html",PTG_transfer=session.get('PTG_transfer', None), session=session)


@app.route("/peg", methods=["POST","GET"])
def peg_generation():
    session['msg'] = None
    session['clr'] = {'sequence': '#FFFFFF', 'edits': '#FFFFFF'}

    if request.method == "POST":
        session['PEG_sequence'] = request.form["sequence"]
        session['PEG_edits'] = request.form["edits"]
        PEG_edits = session['PEG_edits'].split('|')
        PEG_mode = request.form["mode"]

        pegs_list = []
        for edt in PEG_edits:
            pegs_list.append(edt.split(';'))

        edits_list = pegbldr(session['PEG_sequence'], pegs_list, PEG_mode)

        session['PTG_transfer'] = ''
        for c,peg in enumerate(edits_list):
            session['PTG_transfer'] += str(peg[1]) + ';' +  str(peg[2])
            if c < len(edits_list)-1:
                session['PTG_transfer'] += '|'

        return redirect(url_for('sequence'))
    else:
        return render_template("peg_generation.html", PEG_transfer=session.get('PEG_sequence', None), session=session)
      
        
@app.route("/primer_list", methods=["POST","GET"])
def serve_primers():
    csv = 'List of oligos:\n'
    csv += 'oligo_ID,sequence\n'
    positions = len(session['oligo_index'])
    collapsed_index = int(session['oligo_index'])
    oligo_ids = []
    
    for primer in flattn([[i.primer_forward, i.primer_reverse] for i in session['plcstrn'].parts]):
        oligo_ids.append(session['oligo_prefix']+format(collapsed_index,'0'+str(positions)))
        csv += session['oligo_prefix']+format(collapsed_index,'0'+str(positions))+','+primer+'\n'
        collapsed_index += 1
    
    csv += '\nTable of fragments:\n'
    csv += 'fragment_id,fragment_type,forward_primer,Tm_forw,reverse_primer,Tm_rev\n'
    for c,fragment in enumerate(session['plcstrn'].parts):
        csv += session['PTG_name']+'_f'+str(c)+','+fragment.type+','+oligo_ids[c*2]+','+str(np.round(fragment.primer_forward_tm,1))+','+oligo_ids[c*2+1]+','+str(np.round(fragment.primer_reverse_tm, 1))+'\n'
    
    sr = SeqRecord(seq=Seq(session['plcstrn'].sequence, alphabet=IUPAC.ambiguous_dna), name=session['PTG_name'], annotations={'date': date.today().strftime("%d-%b-%Y").upper(), 'topology': 'linear'})
    for ftr in session['plcstrn'].features:
        sr.features.append(ftr)
    gb_json = {}
    gb_json['polycistron'] = polyToJson(session['plcstrn'])
    gb_json['msg'] = session['msg']
    
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
