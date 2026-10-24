from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from urllib.parse import parse_qs
import base64,re,json,datetime
ROOT=Path(__file__).resolve().parent/"browser"
class Handler(BaseHTTPRequestHandler):
 def do_GET(self):
  body=b'<html><body><h1>Local trial evidence saver</h1><form method="post"><label>Filename <input name="filename" aria-label="Filename"></label><label>Base64 content <textarea name="content" aria-label="Base64 content"></textarea></label><button type="submit">Save evidence</button></form></body></html>'
  self.send_response(200);self.send_header("Content-Type","text/html");self.end_headers();self.wfile.write(body)
 def do_POST(self):
  n=int(self.headers.get("Content-Length","0"))
  if n>10000000:self.send_error(413);return
  q=parse_qs(self.rfile.read(n).decode());name=q.get("filename",[""])[0]
  if not re.fullmatch(r"[a-z0-9_-]+\.(png|jpg|txt|json)",name):self.send_error(400);return
  try:data=base64.b64decode(q["content"][0],validate=True)
  except Exception:self.send_error(400);return
  dst=ROOT/name
  if dst.exists():self.send_error(409);return
  if name.endswith(".png") and not data.startswith(b"\x89PNG\r\n\x1a\n"):self.send_error(400);return
  dst.write_bytes(data)
  with (ROOT/"save-log.jsonl").open("a") as f:f.write(json.dumps({"file":name,"bytes":len(data),"utc":datetime.datetime.now(datetime.timezone.utc).isoformat()})+"\n")
  body=("<html><body><h1>Saved "+name+"</h1><p>"+str(len(data))+" bytes</p><a href='/'>Save another</a></body></html>").encode()
  self.send_response(200);self.send_header("Content-Type","text/html");self.end_headers();self.wfile.write(body)
HTTPServer(("127.0.0.1",8767),Handler).serve_forever()
