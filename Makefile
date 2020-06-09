build-docker:
	docker build . -t czbiohub/sc2-msspe:latest

benchmark:
	nextflow run main.nf -profile docker,benchmark
