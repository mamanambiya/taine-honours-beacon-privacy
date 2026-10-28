import os, subprocess, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import real_beacon_attack as rba
ROOT = os.path.expanduser(os.environ.get("POP65_BUILD", "~/pop65_build"))
if not os.path.isdir(ROOT + "/CEU"):
    sys.exit(f"POP65_BUILD not found or not a build directory: {ROOT} (set POP65_BUILD to the folder holding CEU, ESN, ...)")
OUT = os.path.expanduser(os.environ.get("POP65_ATTACK_OUT", ROOT + "/attack"))
os.makedirs(OUT, exist_ok=True)
GROUPS = {"AFR": "ESN GWD LWK MSL YRI AFR".split(), "EUR": "CEU FIN GBR IBS TSI EUR".split()}
CONT = {p: g for g, ps in GROUPS.items() for p in ps}
rd = lambda f: [l.strip() for l in open(f) if l.strip()]
mem = {p: rd(f"{ROOT}/{p}/members.txt") for p in CONT}
ctl = {p: rd(f"{ROOT}/{p}/controls.txt") for p in CONT}
hetdirs = [f"{ROOT}/{p}/het" for p in CONT]
for p in CONT:
    ds = f"{p}_M65_S1"; sites = f"{OUT}/sites_{p}.csv"
    if not os.path.exists(sites) or os.path.getsize(sites) == 0:
        cmd = ["docker", "exec", "beacon-dev-mongodb", "mongoexport", "--db", "beacon_db", "--collection", "variants",
               "--query", json.dumps({"dataset_ids": [ds]}), "--type=csv", "--fields", "reference_name,start,reference_bases,alternate_bases"]
        with open(sites, "w") as f: subprocess.run(cmd, stdout=f, stderr=subprocess.DEVNULL, check=True)
    members = set(mem[p])
    same_pop_controls = sorted(set(ctl[p]) - members)
    continent = set()
    for q in GROUPS[CONT[p]]: continent |= set(mem[q]) | set(ctl[q])
    continent_outsiders = sorted(continent - members)
    for tag, outs in (("controls", same_pop_controls), ("continent", continent_outsiders)):
        mf, of = f"{OUT}/m_{p}.txt", f"{OUT}/o_{p}_{tag}.txt"
        open(mf, "w").write("\n".join(sorted(members)) + "\n"); open(of, "w").write("\n".join(outs) + "\n")
        args = ["--sites", sites, "--dataset", ds, "--members", mf, "--outsiders", of, "--out", f"{OUT}/attack_{p}_{tag}.json"]
        for d in hetdirs: args += ["--het-dir", d]
        t = time.time(); rba._main(args); print(f"  [{p} {tag}] {len(outs)} outsiders, {time.time()-t:.0f}s", flush=True)
print("ATTACK_ALL_DONE", flush=True)
