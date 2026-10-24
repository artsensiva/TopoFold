"""CI helper: fail if a JUnit report contains skipped tests other than an allowed list (exact test names or prefixes).

    python tests/check_junit_skips.py report.xml ALLOWED_PREFIX [ALLOWED_PREFIX ...]
"""
import sys
import xml.etree.ElementTree as ET

report, allowed = sys.argv[1], sys.argv[2:]
cases = list(ET.parse(report).getroot().iter("testcase"))
skipped = [c.get("name") for c in cases if c.find("skipped") is not None]
bad = [s for s in skipped if not any(s.startswith(a) for a in allowed)]
print(f"{len(cases)} test cases, {len(skipped)} skipped, {len(bad)} unexpected skips")
for s in bad:
    print("unexpected skip:", s)
sys.exit(1 if bad or not cases else 0)
