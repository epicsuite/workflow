import yaml

workflow = None
with open("workflow.yaml") as stream:
    try:
        workflow = yaml.load(stream, Loader = yaml.CLoader)
    except yaml.YAMLError as exc:
        print(exc)

ID = 0
for d in workflow['datasets']:
    print('dataset' + str(ID))
    tstep = 0
    for t in d['fastq']:
        print('  timestep' + str(tstep) + ": " + t)
        tstep += 1
    ID += 1

