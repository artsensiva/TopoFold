import os,sys,time,json,subprocess,datetime,pathlib,resource
BASE=pathlib.Path(__file__).resolve().parent
label=sys.argv[1]; command=sys.argv[2:]
env=os.environ.copy()
changes={"VIRTUAL_ENV":str(BASE/".venv"),"PATH":str(BASE/".venv"/"bin")+os.pathsep+env.get("PATH",""),"PIP_CACHE_DIR":str(BASE/"cache"/"pip"),"CARGO_HOME":str(BASE/"cargo"),"CARGO_TARGET_DIR":str(BASE/"target"),"TMPDIR":str(BASE/"tmp"),"MPLCONFIGDIR":str(BASE/"cache"/"matplotlib"),"XDG_CACHE_HOME":str(BASE/"cache"),"PYTHONDONTWRITEBYTECODE":"1","PYTHONPATH":str(BASE/"checkout"/"python"),"PIP_DISABLE_PIP_VERSION_CHECK":"1"}
env.update(changes)
log=BASE/"evidence"/(label+".log")
record={"command":command,"cwd":str(BASE/"checkout"),"environment_overrides":changes,"start_utc":datetime.datetime.now(datetime.timezone.utc).isoformat()}
start=time.monotonic(); usage=resource.getrusage(resource.RUSAGE_CHILDREN)
with log.open("w") as f:
 print("Command:",command,flush=True)
 try:
  p=subprocess.run(command,cwd=BASE/"checkout",env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1200)
  code=p.returncode
 except subprocess.TimeoutExpired:
  code=124; f.write("\nTRIAL COMMAND TIMEOUT after 1200s\n")
record.update(exit_code=code,wall_seconds=time.monotonic()-start,end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
after=resource.getrusage(resource.RUSAGE_CHILDREN)
record.update(cpu_user_seconds=after.ru_utime-usage.ru_utime,cpu_system_seconds=after.ru_stime-usage.ru_stime,max_rss_kib=after.ru_maxrss)
(BASE/"evidence"/(label+".json")).write_text(json.dumps(record,indent=2))
print(json.dumps(record,indent=2)); print(log.read_text()[-18000:]);sys.exit(code)
