docker run --rm  -p 80:80 -p 8889:8888 -p 5000:5000 -e IP_MAC=$(ipconfig getifaddr en0) -v $(pwd)/polygen_scripts/:/scripts/ -w /scripts/ -it uurquiza/polygen:latest
