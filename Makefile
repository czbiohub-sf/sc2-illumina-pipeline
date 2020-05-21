build-docker:
	docker build . -t czbiohub/sc2-msspe:latest

benchmark:
	nextflow run call_consensus.nf -profile docker,benchmark
