import os, sys, random, json, urllib.request, urllib.parse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import real_beacon_attack as rba
R = os.path.expanduser(os.environ.get("POP65_BUILD", "~/pop65_build"))
rng = random.Random(1)
tot = bad = 0
for pop, people in (("LWK", ["NA19445", None]), ("EUR", [None, None])):
    sites = rba.load_beacon_sites(f"{os.path.expanduser(os.environ.get('POP65_ATTACK_OUT', R + '/attack'))}/sites_{pop}.csv")
    members = [l.strip() for l in open(f"{R}/{pop}/members.txt")]
    controls = [l.strip() for l in open(f"{R}/{pop}/controls.txt")]
    # one outsider (a control) and one member per beacon
    for who, label in ((controls[0] if people[0] is None else people[0], "outsider(control)"), (members[0], "member")):
        path = None
        for d in os.listdir(R):
            p = f"{R}/{d}/het/{who}.tsv"
            if os.path.exists(p): path = p; break
        rows = [l.rstrip("\n").split("\t") for l in open(path)]
        for chrom, pos, ref, alt in rng.sample(rows, 100):
            q = urllib.parse.urlencode({"assemblyId": "GRCh38", "referenceName": chrom.replace("chr", ""), "start": int(pos) - 1,
                                        "referenceBases": ref, "alternateBases": alt, "datasetIds": f"{pop}_M65_S1"})
            r = json.load(urllib.request.urlopen("http://127.0.0.1:8086/api/g_variants?" + q))
            api = bool(r["responseSummary"]["exists"]); local = rba.site_key(chrom, pos, ref, alt) in sites
            tot += 1; bad += (api != local)
        print(f"{pop} {label} {who}: 100 queries done")
print("API vs set lookup: queries", tot, "disagreements", bad)
