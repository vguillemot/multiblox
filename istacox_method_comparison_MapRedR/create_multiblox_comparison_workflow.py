import sys
import json
from soma_workflow.client import Job, Workflow, Helper
# args
all_jobs_file = sys.argv[1]
method_file = sys.argv[2]
workflow_file = sys.argv[3]

with open(method_file, 'r') as f:
    method_cfg = json.load(f)
try:
    test = method_cfg["model.inner.reducer"]
except:
    method_cfg["model.inner.reducer"] = "None"

jobs = []
red_jobs = {}
dep_dic = {}
# read file
with open(all_jobs_file, 'r') as f:
    for line in f:
        cmd = line.split(" ")[:-1]

        if(cmd[1].endswith("multiblox_comparison.mapper.R")):
            print "multiblox_comparison.mapper"
            design = cmd[-2] #cmd[2].split('_')[-2]
            outer = cmd[6]
            key = outer+design
            if(not dep_dic.has_key(key)):
                dep_dic[key] = []                
            cur_job = Job(command=cmd, name="Map %s outer %s inner %s" %
                          (design, outer, cmd[5]))
            dep_dic[key].append(cur_job)
        elif(cmd[1].endswith("multiblox_comparison.inner.reducer.R")):
            print "multiblox_comparison.inner.reducer"
            design = cmd[-1] #.split('_')[-2]
            outer = cmd[5]
            key = outer+design
            if(not dep_dic.has_key(key)):
                dep_dic[key] = []
            cur_job = Job(command=cmd, name="Reduce %s outer %s" %
                          (design, outer))
            red_jobs[key] = cur_job
        elif(cmd[1].endswith("multiblox_comparison.outer.reducer.R")):
            print "multiblox_comparison.outer.reducer"
            cur_job = Job(command=cmd, name="Final Reduce")
            final_red = cur_job
        elif(cmd[1].endswith("coxnet.mapper.R")):
            print "coxnet mapper"
            cur_job = Job(command=cmd, name="Coxnet")
            glmnet_job = cur_job
            pass
        else:
            raise Exception("Unknown task, abort...")
        jobs.append(cur_job)            
dependencies = []
for k, j in red_jobs.items():
    parent_list = dep_dic[k]
    for p in parent_list:
        dependencies.append((p, j))
    dependencies.append((j, final_red))
dependencies.append((glmnet_job, final_red))
workflow = Workflow(jobs=jobs,
                    dependencies=dependencies)

# save the workflow into a file
Helper.serialize(workflow_file, workflow)
