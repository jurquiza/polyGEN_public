from flask import Flask

application = Flask(__name__)


#we need to genearate the correct credentials for logging into synbiohub.
@application.route('/')
def hello_world():
    return 'hello, world I am polygen!!'

if __name__=='__main__':
    application.run(debug=True, host='0.0.0.0', port=5000) ##host 0.0.0.0 is important for docker container
