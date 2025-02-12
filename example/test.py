import yaml

workflow = None
with open("workflow.yaml") as stream:
    try:
        workflow = yaml.load(stream, Loader = yaml.CLoader)
    except yaml.YAMLError as exc:
        print(exc)

print(workflow['datasets'][0])
print(workflow['datasets'][0]['fastq'][0])

