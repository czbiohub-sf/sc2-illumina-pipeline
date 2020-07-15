build-docker:
	docker build . -t czbiohub/sc2-msspe:latest

benchmark:
	nextflow run czbiohub/sc2-illumina-pipeline --profile docker,benchmark
