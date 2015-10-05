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

        if(cmd[1].endswith(method_cfg["model.mapper"])):
            outer = cmd[6] if method_cfg["model.selection"] else cmd[4]
            key = outer
            if(not dep_dic.has_key(key)):
                dep_dic[key] = []
            if method_cfg["model.selection"]:
                jname = "Map outer %s inner %s" % (outer, cmd[5])  
            else:
                jname = "Map outer %s" % outer
            cur_job = Job(command=cmd, name=jname)
            dep_dic[key].append(cur_job)
        elif(cmd[1].endswith(method_cfg["model.inner.reducer"])):
            outer = cmd[5]
            key = outer
            if(not dep_dic.has_key(key)):
                dep_dic[key] = []
            cur_job = Job(command=cmd, name="Reduce inner %s" %
                          (outer))
            red_jobs[key] = cur_job
        elif(cmd[1].endswith(method_cfg["model.outer.reducer"])):
            cur_job = Job(command=cmd, name="Final Reduce")
            final_red = cur_job
        else:
            raise Exception("Unknown task, abort...")
        jobs.append(cur_job)

dependencies = []
if (method_cfg["model.selection"]):
    for k, j in red_jobs.items():
        parent_list = dep_dic[k]
        for p in parent_list:
            dependencies.append((p, j))
        dependencies.append((j, final_red))
else:
   for _, j in dep_dic.items():
       dependencies.append((j[0], final_red))

workflow = Workflow(jobs=jobs,
                    dependencies=dependencies)

# save the workflow into a file
Helper.serialize(workflow_file, workflow)
