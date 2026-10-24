from pathlib import Path
import os,sys,json,time,subprocess,datetime,resource
b=Path(__file__).resolve().parent; label=sys.argv[1]; guide=sys.argv[2]; command=sys.argv[3:]
e=b/'evidence'; env=os.environ.copy(); overrides={'CARGO_HOME':str(b/'cargo'),'CARGO_TARGET_DIR':str(b/'target'),'PIP_CACHE_DIR':str(b/'cache/pip'),'XDG_CACHE_HOME':str(b/'cache'),'MPLCONFIGDIR':str(b/'cache/matplotlib'),'TMPDIR':str(b/'tmp'),'PYTHONDONTWRITEBYTECODE':'1','PIP_DISABLE_PIP_VERSION_CHECK':'1'};env.update(overrides)
rec={'guide_command':guide,'executed_argv':command,'cwd':str(b/'checkout'),'environment_overrides':overrides,'start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
t=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
with (e/(label+'.log')).open('w') as f:
 try:p=subprocess.run(command,cwd=b/'checkout',env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1200);code=p.returncode
 except subprocess.TimeoutExpired:code=124;f.write('TIMEOUT after 1200 s\n')
after=resource.getrusage(resource.RUSAGE_CHILDREN)
rec.update(exit_code=code,end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),wall_seconds=time.monotonic()-t,cpu_user_seconds=after.ru_utime-before.ru_utime,cpu_system_seconds=after.ru_stime-before.ru_stime,peak_rss_kib=after.ru_maxrss)
(e/(label+'.json')).write_text(json.dumps(rec,indent=2)+'\n'); print(json.dumps(rec,indent=2)); print((e/(label+'.log')).read_text()[-10000:]);sys.exit(code)
