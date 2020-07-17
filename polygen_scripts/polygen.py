from flask import Flask, redirect, url_for, render_template, request

app = Flask(__name__)

@app.route("/")
def home():
    return render_template("index.html", content="Testing")

@app.route("/sequence", methods=["POST","GET"])
def sequence():
    if request.method == "POST":
        sequence = request.form["sequence"]
        return redirect(url_for("success"))
    else:
        return render_template("sequence.html")

@app.route("/sucess")
def success():
    return "Sucess"


if __name__=='__main__':
    app.run(debug=True, host='0.0.0.0', port=5000) ##host 0.0.0.0 is important for docker container
