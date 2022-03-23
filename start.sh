docker run --rm -p 8889:8888 -p 5000:5000 -v $(pwd)/polygen_scripts/:/scripts/ -w /scripts/ -it polygen:test
