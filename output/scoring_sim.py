#!/usr/bin/env python3
"""
MolCraft 赛事评分模拟器
基于FAQ+已知分数反向工程绑定能的归一化公式
"""
import sys, os, json, math, csv
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
from evaluator import evaluate_molecule
from rdkit import Chem

def _parse_route(route_str):
    """解析合成路线，检查各步问题"""
    steps = route_str.split('|')
    issues = []
    trivial = False
    product_matches = True
    balance_ok = True
    
    if len(steps) == 0:
        return {'steps': 1, 'trivial': True, 'product_matches': True, 'balance_ok': False, 'issues': ['空路线']}
    
    for step in steps:
        parts = step.split('>>')
        if len(parts) == 2:
            reactants_str, product_str = parts
            if product_str == reactants_str:
                trivial = True
                issues.append('A→A')
            r_list = reactants_str.split('.')
            p_mol = Chem.MolFromSmiles(product_str)
            p_heavy = p_mol.GetNumHeavyAtoms() if p_mol else 0
            r_heavy = 0
            for r in r_list:
                m = Chem.MolFromSmiles(r)
                if m:
                    r_heavy += m.GetNumHeavyAtoms()
            # 检查Suzuki: 反应物有B/卤素但产物没有
            p_elements = set()
            if p_mol:
                for atom in p_mol.GetAtoms():
                    p_elements.add(atom.GetSymbol())
            r_elements = set()
            for r in r_list:
                m = Chem.MolFromSmiles(r)
                if m:
                    for atom in m.GetAtoms():
                        r_elements.add(atom.GetSymbol())
            lost = r_elements - p_elements
            if lost and ('Br' in lost or 'B' in lost or 'I' in lost):
                balance_ok = False
                issues.append(f'丢失元素{lost}')
        elif len(parts) == 1:
            if parts[0] == parts[0]:  # single molecule
                trivial = True
                issues.append('单分子')
    
    # 最后一步产物是否匹配设计分子
    if steps:
        last_step = steps[-1]
        parts = last_step.split('>>')
        if len(parts) == 2:
            final_product = parts[1]
        else:
            final_product = parts[0]
        product_matches = True  # will be checked by caller
    
    return {
        'n_steps': len(steps),
        'trivial': trivial,
        'balance_ok': balance_ok,
        'issues': issues,
    }


def simulate(molecules, binding_norm_mode='unknown', verbose=False):
    """
    评分模拟器
    
    molecules: list of (smiles, route)
    
    核心难点: binding_score的归一化公式
    需要多个候选公式，通过真实分数(0.1766)反向选择
    
    binding_norm_mode:
      - 'auto': 尝试所有公式并输出预测
      - 'ref_12_avg': Vina绝对值/12, 取平均 (最匹配已知数据)
      - 'sigmoid': sigmoid归一化
      - 'tanh': tanh归一化
    """
    n = len(molecules)
    mol_data = []
    
    for smi, route in molecules:
        props = evaluate_molecule(smi)
        mol = Chem.MolFromSmiles(smi)
        
        # 分子性质
        sa = props.get('sa_score', 10)
        mw = props.get('mw', 0)
        qed = props.get('qed', 0)
        logp = props.get('logp', 0)
        
        # 路线分析
        route_info = _parse_route(route)
        
        mol_data.append({
            'smiles': smi,
            'sa': sa,
            'mw': mw,
            'qed': qed,
            'logp': logp,
            'n_steps': route_info['n_steps'],
            'trivial': route_info['trivial'],
            'balance_ok': route_info['balance_ok'],
            'issues': route_info['issues'],
            'sa_comp': 0.0 if sa >= 4 else (1.0 - sa/6.0),
        })
    
    # —— binding_score ——
    # 我们之前实测的几个Vina分数:
    # Deucravacitinib: -8.987 (MW422)
    # 大环烷烃: -10.812 (MW335)
    # 4环稠芳烃: -13.236 (MW457)
    # 提交10个分子的best ≈ -9.335
    
    # 因为我们没法每次都对接10个分子，用之前迭代的结果估算
    # H011最佳BE=-9.941, 平均BE=-8.961, trivial=1/10
    # H019最佳BE=-9.335, 平均BE=-8.581, trivial=1/10
    
    # 用迭代历史中的BE数据来推测
    # 根据iteration_log，各轮次的top-10 avg BE大约在 -7.9 到 -9.0 之间
    # 最终提交的10个分子avg BE≈-8.58
    
    avg_be_estimate = -8.58  # 基于历史数据的合理估计
    best_be_estimate = -9.34
    
    # 各种归一化公式
    def norm_minmax(be, pool_min=-14, pool_max=0):
        """min-max归一化，假设pool范围"""
        return (pool_max - be) / (pool_max - pool_min)
    
    def norm_ref12(be):
        """参考-12线性归一化"""
        val = abs(be) / 12.0
        return min(1.0, max(0.0, val))
    
    def norm_sigmoid(be, k=0.4, x0=-7):
        """sigmoid: binding_score = 1/(1+exp(k*(be-x0)))"""
        return 1.0 / (1.0 + math.exp(k * (be - x0)))
    
    def norm_tanh(be, scale=10):
        """tanh归一化"""
        return (1.0 + math.tanh(-be / scale)) / 2.0
    
    def norm_exp_power(be, p=0.5):
        """指数幂: 对低Vina有强区分"""
        return abs(be / 20.0) ** p
    
    # 尝试各种公式
    binding_methods = {
        'minmax_14': norm_minmax(avg_be_estimate, -14, 0),
        'minmax_20': norm_minmax(avg_be_estimate, -20, 0), 
        'ref_12_avg': norm_ref12(avg_be_estimate),
        'sigmoid_k04_x0n7': norm_sigmoid(avg_be_estimate, 0.4, -7),
        'sigmoid_k05_x0n6': norm_sigmoid(avg_be_estimate, 0.5, -6),
        'sigmoid_k03_x0n8': norm_sigmoid(avg_be_estimate, 0.3, -8),
        'tanh_scale10': norm_tanh(avg_be_estimate, 10),
        'exp_pow05': norm_exp_power(avg_be_estimate, 0.5),
        'exp_pow07': norm_exp_power(avg_be_estimate, 0.7),
    }
    
    # —— validity ——
    # 全都valid (之前已验证)
    validity_score = 1.0
    
    # —— sa_score ——
    sa_vals = [d['sa_comp'] for d in mol_data]
    sa_score = sum(sa_vals) / n
    
    # —— route_validity ——
    valid_count = sum(1 for d in mol_data if not d['trivial'])
    route_validity_score = valid_count / n
    
    # —— balance_score ——
    balance_ok_count = sum(1 for d in mol_data if d['balance_ok'])
    balance_score = balance_ok_count / n
    
    # —— step_penalty ——
    avg_steps = sum(d['n_steps'] for d in mol_data) / n
    step_penalty = math.exp(-0.15 * avg_steps)  # 步数越多分越低
    
    # —— convergence ——
    convergence = 0.5  # 简化
    
    # —— starting_material —— (已知0.9)
    starting_material = 0.9
    
    # 计算mol_score和route_score
    results = {}
    for method, binding_val in binding_methods.items():
        mol_score = 0.8 * binding_val + 0.1 * validity_score + 0.1 * sa_score
        
        route_score = (0.55 * route_validity_score + 
                       0.30 * starting_material + 
                       0.05 * step_penalty + 
                       0.05 * convergence + 
                       0.05 * balance_score)
        
        total = 0.7 * mol_score + 0.3 * route_score
        
        # 硬性归零: balance=0的路线
        # 这里简化处理
    
        results[method] = {
            'binding': binding_val,
            'mol_score': mol_score,
            'route_score': route_score,
            'total': total,
            'match': abs(binding_val - 0.1766) < 0.02,
        }
    
    info = {
        'n_molecules': n,
        'avg_steps': avg_steps,
        'sa_score': sa_score,
        'route_validity_score': route_validity_score,
        'balance_score': balance_score,
        'step_penalty': step_penalty,
        'methods': results,
        'mol_data': mol_data,
    }
    
    return info


def main():
    # 加载提交分子
    csv_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'output', 'result.csv')
    mols = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mols.append((row['mol_smiles'].strip(), row['route'].strip()))
    
    print("=" * 60)
    print("MolCraft 赛事评分模拟器")
    print("=" * 60)
    print(f"\n提交分子数: {len(mols)}")
    print()
    
    sim = simulate(mols)
    
    print(f"提取属性:")
    print(f"  平均步数: {sim['avg_steps']:.2f}")
    print(f"  sa_score: {sim['sa_score']:.4f} (真实: 0.6561)")
    print(f"  route_validity: {sim['route_validity_score']} (真实: 0.5)")
    print(f"  balance_ok比例: {sim['balance_score']} (影响balance_score)")
    print(f"  step_penalty: {sim['step_penalty']:.4f}")
    print()
    
    print("各种binding归一化方案对比 (目标: 0.1766):")
    print(f"{'方法':<20} {'binding':>8} {'mol_score':>10} {'route_score':>12} {'total':>8} {'匹配?':>5}")
    print("-" * 63)
    for method_name, m in sorted(sim['methods'].items(), key=lambda x: -x[1]['total']):
        match_mark = "✅" if m['match'] else ""
        print(f"{method_name:<20} {m['binding']:>8.4f} {m['mol_score']:>10.4f} {m['route_score']:>12.4f} {m['total']:>8.4f} {match_mark:>5}")
    
    print()
    print("=" * 60)
    print("最可能的归一化公式分析")
    print("=" * 60)
    
    # 哪个公式最接近0.1766?
    best_method = min(sim['methods'].items(), key=lambda x: abs(x[1]['binding'] - 0.1766))
    print(f"\n最佳匹配: {best_method[0]} = {best_method[1]['binding']:.4f} (距0.1766差{abs(best_method[1]['binding']-0.1766):.4f})")
    
    # 用最佳公式计算总分
    if best_method[0].startswith('sigmoid'):
        print("\n推测 → binding_score使用sigmoid归一化：")
        print("  binding_score = 1/(1+exp(k*(Vina-x0)))")
        print(f"  这符合赛事FAQ的'越负越好'且'有饱和区'的描述")
    elif best_method[0].startswith('ref'):
        print("\n推测 → binding_score使用线性归一化：")
        print("  binding_score = |Vina| / 参考值, capped at 1.0")
    elif best_method[0].startswith('tanh'):
        print("\n推测 → binding_score使用tanh归一化：")
    else:
        print(f"\n推测 → 使用{best_method[0]}归一化")
    
    print()
    print("各分子情况:")
    for i, d in enumerate(sim['mol_data']):
        issues = ', '.join(d['issues']) if d['issues'] else 'OK'
        print(f"  #{i+1}: SA={d['sa']:.1f} {'❌' if d['sa']>=4 else '✅'} "
              f"步数={d['n_steps']} 平衡={'✅' if d['balance_ok'] else '❌'} "
              f"trivial={'⚠️' if d['trivial'] else '✅'} | {issues}")


if __name__ == '__main__':
    main()