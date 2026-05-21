"""
Regression test for the SMARTS mis-match fix.
Tests that:
1. Problem molecules are now correctly rejected
2. Valid molecules still produce valid routes
"""
import sys
sys.path.insert(0, '/home/dministrator/Lab/clones/molcraft-agent/src')
from synthesis_v2 import plan_synthesis_v2, RETRO_RULES

print("=" * 80)
print("REGRESSION TEST: Problem molecules should be rejected/fixed")
print("=" * 80)

# The full target molecules from result.csv
problem_targets = [
    ("#2 Cc1ccc(-c2cncc3ccccc23)cc1", "Cc1ccc(-c2cncc3ccccc23)cc1"),
    ("#3 FCc1ncc2cccc(-c3ccccc3)c2n1", "FCc1ncc2cccc(-c3ccccc3)c2n1"),
    ("#4 c1ccc(-c2cccc3ncncc23)cc1", "c1ccc(-c2cccc3ncncc23)cc1"),
]

for name, smi in problem_targets:
    print(f"\n  {name}: {smi}")
    result = plan_synthesis_v2(smi)
    route = result.get("route", "N/A")
    steps = result.get("steps", 0)
    trivial = result.get("trivial", False)
    
    # Check if route contains the problematic "C=O>>..." pattern in the first step
    has_co_issue = ".C=O>>" in route or ".C=O>>" in route
    
    print(f"    Route: {route[:80]}...")
    print(f"    Steps: {steps}, Trivial: {trivial}")
    print(f"    Has C=O issue: {has_co_issue}")
    
    if has_co_issue:
        print(f"    ❌ STILL HAS PROBLEM: C=O in reactants without matching heteroatoms")
    else:
        print(f"    ✅ No C=O issue (either rejected or different route)")

print("\n" + "=" * 80)
print("REGRESSION TEST: Valid molecules should still work")
print("=" * 80)

valid_targets = [
    ("Simple Suzuki biphenyl", "c1ccc(-c2ccccc2)cc1"),
    ("Benzamide", "O=C(Nc1ccccc1)C"),
    ("Simple quinazoline", "c1ccc2ncncc2c1"),
    ("Simple isoquinoline", "c1ccc2ccnc2c1"),
    ("Benzyl alcohol", "OCc1ccccc1"),
]

for name, smi in valid_targets:
    print(f"\n  {name}: {smi}")
    result = plan_synthesis_v2(smi)
    success = result.get("success", False)
    route = result.get("route", "N/A")
    trivial = result.get("trivial", False)
    steps = result.get("steps", 0)
    
    print(f"    Success: {success}")
    print(f"    Route: {route[:80]}")
    print(f"    Steps: {steps}, Trivial: {trivial}")
    
    if success and not trivial and steps > 0:
        print(f"    ✅ Valid non-trivial route")
    elif trivial:
        print(f"    ⚠️ Trivial route (might need more rules)")
    elif not success:
        print(f"    ❌ Failed to find any route")

print("\n" + "=" * 80)
print("SMARTS RULE VERIFICATION: Updated quinazoline SMARTS")
print("=" * 80)

from rdkit import Chem

# Find the quinazoline rule
for smarts_pat, retro_smarts in RETRO_RULES:
    if "ncnc" in smarts_pat:
        print(f"\n  SMARTS: {smarts_pat}")
        print(f"  Retro: {retro_smarts}")
        
        # Test on true quinazoline
        qz = Chem.MolFromSmiles("c1ccc2ncncc2c1")
        qz_patt = Chem.MolFromSmarts(smarts_pat)
        print(f"  Matches true quinazoline: {qz.HasSubstructMatch(qz_patt)}")
        
        # Test on problem molecule #3 intermediate
        mol3 = Chem.MolFromSmiles("FCc1ncc2cccc(Br)c2n1")
        print(f"  Matches naphthyridine #3: {mol3.HasSubstructMatch(qz_patt)}")

print("\n" + "=" * 80)
print("TEST COMPLETE")
print("=" * 80)