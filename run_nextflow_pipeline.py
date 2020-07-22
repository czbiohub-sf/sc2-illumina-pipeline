
import argparse
import subprocess
import tempfile
import json

import boto3
from botocore.exceptions import ClientError


def get_secret():

    secret_name = "czb-covid-pipeline-weblogs-endpoint"
    region_name = "us-west-2"

    # Create a Secrets Manager client
    session = boto3.session.Session()
    client = session.client(
        service_name='secretsmanager',
        region_name=region_name
    )

    try:
        get_secret_value_response = client.get_secret_value(
            SecretId=secret_name
        )
    except ClientError as e:
        raise e

    secret = json.loads(get_secret_value_response['SecretString'])
    return secret["weblogs-http-endpoint"]


def replace_string(filename, old_string, new_string):
    with open(filename) as f:
        s = f.read()
        if old_string not in s:
            print('"{old_string}" not found in {filename}.'.format(**locals()))
        else:
            s = s.replace(old_string, new_string)
            return s


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run nextflow pipeline on Biohub infra, and in a way that covidhub db can populate from the results."
    )

    args, remaining_nextflow_args = parser.parse_known_args()

    secret = get_secret()
    nextflow_config = replace_string("conf/awsbatch.config", "$WEBLOGS_ENDPOINT", secret)

    cmd = ["nextflow", "run", "czbiohub/sc2-illumina-pipeline", "-resume"]

    with tempfile.NamedTemporaryFile(mode="w") as f:
        print(nextflow_config, file=f)
        cmd += ["-c", f.name]
        cmd += remaining_nextflow_args
        subprocess.run(cmd)
