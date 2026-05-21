#!/usr/bin/env python3
"""Generate and score different route format versions for platform compatibility testing."""
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'src'))

INPUT_CSV = os.path.join(os.path.dirname(__file__), 'result.csv')
SINGLE_STEP_OUT = os.path.join(os.path.dirname(__file__), 'single_step_result.csv')
NO_BYPROD_OUT = os.path.join(os.path.dirname(__file__), 'single_step_nobyproduct_result.csv')

def parse_step(step: str):
    """Parse a reaction step into reactants and products."""
    if '>>' not in step:
        return [], []
    r, p = step.split('>>', 1)
    reactants = [x.strip() for x in r.split('.') if x.strip()]
    products = [x.strip() for x in p.split('.') if x.strip()]
    return reactants, products

def strip_byproducts(route: str) -> str:
    """Remove byproducts from product side, keeping only the first (main) product."""
    reactants, products = parse_step(route)
    if not products:
        return route
    # Keep only the first product (the actual target molecule)
    return f"{'.'.join(reactants)}>>{products[0]}"

def main():
    # Read input
    rows = []
    with open(INPUT_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    # Generate single-step + byproducts version
    single_rows = []
    for row in rows:
        mol_smiles = row.get('mol_smiles', '').strip()
        route = row.get('route', '').strip()
        # Extract last step
        steps = route.split(' | ')
        last_step = steps[-1].strip()
        single_rows.append({'mol_smiles': mol_smiles, 'route': last_step})
    
    with open(SINGLE_STEP_OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['mol_smiles', 'route'])
        w.writeheader()
        w.writerows(single_rows)
    print(f"✅ Single-step (with byproducts): {SINGLE_STEP_OUT}")

    # Generate single-step + NO byproducts version
    nobyprod_rows = []
    for row in rows:
        mol_smiles = row.get('mol_smiles', '').strip()
        route = row.get('route', '').strip()
        steps = route.split(' | ')
        last_step = steps[-1].strip()
        # Remove byproducts from product side
        no_bp_route = strip_byproducts(last_step)
        nobyprod_rows.append({'mol_smiles': mol_smiles, 'route': no_bp_route})
    
    with open(NO_BYPROD_OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['mol_smiles', 'route'])
        w.writeheader()
        w.writerows(nobyprod_rows)
    print(f"✅ Single-step no-byproduct: {NO_BYPROD_OUT}")

    # Run local scoring on all three versions
    print()
    print("=" * 60)
    print("LOCAL SCORING COMPARISON (using scorer.py)")
    print("=" * 60)
    
    from src.scorer import score_csv
    
    for label, path in [("Original (multi-step)", INPUT_CSV),
                         ("Single-step (with byproducts)", SINGLE_STEP_OUT),
                         ("Single-step (no byproducts)", NO_BYPROD_OUT)]:
        try:
            result = score_csv(path)
            print(f"\n{label}:")
            print(f"  route_validity_score: {result['route_validity_score']}")
            print(f"  balance_score (avg):  {result.get('route_score'):.6f} (route_score component includes balance)")
            print(f"  step_penalty_score:  {result['starting_material_availability_score']:.6f}")
            route_validity = result.get('route_validity_score', 0)
            n_valid = int(route_validity * 10)
            print(f"  → {n_valid}/10 routes pass platform validity check (predicted)")
        except Exception as e:
            print(f"\n{label}: ERROR - {e}")
    print()

    # Show the no-byproduct routes
    print("=" * 60)
    print("NO-BYPRODUCT ROUTES PREVIEW")
    print("=" * 60)
    for i, row in enumerate(nobyprod_rows):
        print(f"  #{i+1}: {row['route'][:100]}...")

if __name__ == '__main__':
    main()