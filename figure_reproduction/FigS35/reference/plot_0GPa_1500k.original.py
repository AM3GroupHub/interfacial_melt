# -*- coding: utf-8 -*-
"""
能量-组分关系脚本（端元归一化 + 5点三阶拟合 + 画拟合曲线 + 计算中点曲率并标注）
修改：曲线颜色根据曲率正负分段变化（负红正蓝）
"""

import re
import numpy as np
import matplotlib.pyplot as plt
from math import gcd
from functools import reduce

# ==================== 参数设置 ====================
FILES = [
    'Fe324B36.txt',
    'Fe128B32.txt',
    'Fe252B108.txt',
    'Fe216B144.txt',
    'Fe200B200.txt',
    'Fe160B240.txt',
    'Fe135B405.txt',
    'Fe120B280.txt',
    'Fe81B459.txt',
    'Fe72B288.txt',
    'Fe54B486.txt'     
]
REFERENCE_FORMULA = 'B'
TEMPERATURE = 1500
OUTPUT_PREFIX = 'energy_composition'

POLY_DEGREE = 3          # 固定三阶
N_FIT_SAMPLES = 300      # 拟合曲线采样点数
KAPPA_FMT = ".3e"        # 曲率标注格式
# =================================================


def parse_formula(formula):
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)
    composition = {}
    for element, count in matches:
        count = int(count) if count else 1
        composition[element] = composition.get(element, 0) + count
    return composition


def read_energy_file(filename):
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    temp = float(parts[0])
                    pot_eng = float(parts[1]) if parts[1] != 'N/A' else np.nan
                    enthalpy = float(parts[2]) if parts[2] != 'N/A' else np.nan
                    data[temp] = (pot_eng, enthalpy)
                except ValueError:
                    continue
    return data


def gcd_of_list(nums):
    nums = [int(n) for n in nums if int(n) != 0]
    return reduce(gcd, nums) if nums else 1


def composition_to_formula(comp):
    parts = []
    for el in sorted(comp.keys()):
        n = int(comp[el])
        parts.append(f"{el}{'' if n == 1 else n}")
    return "".join(parts)


def reduce_formula_with_factor(compound_formula):
    comp = parse_formula(compound_formula)
    g = gcd_of_list(list(comp.values()))
    comp_min = {k: v // g for k, v in comp.items()}
    return comp_min, g, composition_to_formula(comp_min)


def endmember_decompose_and_divisor(compound_formula, reference_formula):
    """
    最简式 = n_ref * reference + n_rem * remainder_min
    x_ref = n_ref/(n_ref+n_rem)
    divisor = g*(n_ref+n_rem)
    """
    comp_min, g, comp_min_str = reduce_formula_with_factor(compound_formula)
    ref = parse_formula(reference_formula)

    n_ref = float('inf')
    for el, cnt in ref.items():
        if el not in comp_min:
            n_ref = 0
            break
        n_ref = min(n_ref, comp_min[el] // cnt)
    n_ref = int(n_ref if n_ref != float('inf') else 0)

    rem = comp_min.copy()
    for el, cnt in ref.items():
        rem[el] = rem.get(el, 0) - n_ref * cnt
        if rem.get(el, 0) == 0:
            rem.pop(el, None)
    rem = {k: v for k, v in rem.items() if v > 0}

    if not rem:
        n_rem = 0
        remainder_str = ""
    else:
        n_rem = gcd_of_list(list(rem.values()))
        rem_min = {k: v // n_rem for k, v in rem.items()}
        remainder_str = composition_to_formula(rem_min)

    total = n_ref + n_rem
    x_ref = n_ref / total if total > 0 else 0.0
    divisor = g * total if total > 0 else np.nan

    return x_ref, divisor, remainder_str, n_ref, n_rem, g, comp_min_str


def cubic_fit_and_mid_curvature(x, y):
    """
    输入必须是 5 个点，用三阶多项式拟合 y(x)，并计算中间点(第3个)曲率：
      κ = y''(x0) / (1 + (y'(x0))^2)^(3/2)
    返回：
      p(多项式), kappa_mid, x0, y0
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) != 5 or len(y) != 5:
        raise ValueError("该算法要求输入恰好 5 个点")

    coeffs = np.polyfit(x, y, deg=3)
    p = np.poly1d(coeffs)
    dp = np.polyder(p, 1)
    ddp = np.polyder(p, 2)

    mid = 2
    x0 = float(x[mid])
    y0 = float(y[mid])

    y1 = float(dp(x0))
    y2 = float(ddp(x0))

    kappa = y2 / ((1.0 + y1**2) ** 1.5)
    return p, kappa, x0, y0


def annotate_single_kappa(ax, x0, y0, kappa, fmt):
    y0min, y0max = ax.get_ylim()
    yspan = y0max - y0min
    dy = 0.03 * yspan if np.isfinite(yspan) and yspan != 0 else 0.1
    
    col = "tab:red" if kappa > 0 else ("tab:blue" if kappa < 0 else "black")
    
    ax.text(
        x0, y0 + dy, format(kappa, fmt),
        fontsize=10, color=col, ha="center", va="bottom",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.6, linewidth=0.6)
    )


def plot_segmented_curve_by_curvature(ax, p, x_range, n_samples=N_FIT_SAMPLES):
    """
    将拟合曲线按曲率正负分段绘制
    p: 多项式对象
    x_range: (xmin, xmax)
    """
    xs = np.linspace(x_range[0], x_range[1], n_samples)
    ys = p(xs)
    
    # 计算二阶导数（曲率的分子部分）
    dp = np.polyder(p, 1)
    ddp = np.polyder(p, 2)
    
    # 计算每个点的曲率
    y_prime = dp(xs)
    y_double_prime = ddp(xs)
    curvatures = y_double_prime / ((1.0 + y_prime**2) ** 1.5)
    
    # 找到曲率符号变化的点
    sign_changes = np.where(np.diff(np.sign(curvatures)))[0]
    
    # 分段绘制
    segments = []
    start_idx = 0
    
    for change_idx in sign_changes:
        segments.append((start_idx, change_idx + 1))
        start_idx = change_idx + 1
    segments.append((start_idx, len(xs)))
    
    # 绘制每一段
    for start, end in segments:
        if start >= end:
            continue
        
        x_seg = xs[start:end]
        y_seg = ys[start:end]
        
        # 使用该段中点的曲率判断颜色
        mid_idx = (start + end) // 2
        if mid_idx >= len(curvatures):
            mid_idx = len(curvatures) - 1
        
        kappa_mid = curvatures[mid_idx]
        
        if kappa_mid < 0:
            color = 'red'
        elif kappa_mid > 0:
            color = 'blue'
        else:
            color = 'gray'
        
        ax.plot(x_seg, y_seg, color=color, linewidth=2.5, alpha=0.9)
    
    # 添加图例（只添加一次）
    from matplotlib.lines import Line2D
    custom_lines = [
        Line2D([0], [0], color='blue', lw=2.5),
        Line2D([0], [0], color='red', lw=2.5)
    ]
    ax.legend(custom_lines, ['κ > 0 (Convex)', 'κ < 0 (Concave)'], fontsize=10, frameon=True)



def plot_one_panel(ax, x, y, color_points, color_line_default, title, ylabel, reference_formula):
    """
    绘制单个面板，包含分段颜色的拟合曲线（支持4个及以上数据点，不再计算和标注中点曲率）
    """
    ax.scatter(x, y, s=90, alpha=0.80, edgecolors='black', linewidth=1.0, color=color_points)
    ax.plot(x, y, '-', alpha=0.25, color=color_points)

    # 三阶多项式拟合至少需要4个点
    if len(x) >= 4:
        # 直接进行三阶多项式拟合获取多项式对象 p
        coeffs = np.polyfit(x, y, 6)
        p = np.poly1d(coeffs)

        # 使用分段绘制函数绘制拟合曲线
        plot_segmented_curve_by_curvature(ax, p, (np.min(x), np.max(x)))

        print(f"{title}：已完成三阶拟合并绘制曲线")
    else:
        print(f"{title}：点数少于4个，数据不足，未进行三阶拟合")

    ax.set_xlabel(f'Mole ratio of {reference_formula} (x)', fontsize=12, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)


def plot_energy_composition(files, reference_formula, temperature, output_prefix):
    data_pot = []  # (x_ref, pot_norm)
    data_ent = []  # (x_ref, ent_norm)

    P0 = -8.043070824
    P1 = -6.434663479166667
    
    ent0 = -7.850557452
    ent1 = -6.2431970833333335

    for filename in files:
        compound = filename.replace('.txt', '')
        energy_data = read_energy_file(filename)
        if temperature not in energy_data:
            print(f"警告: {filename} 中没有温度 {temperature} K 的数据")
            continue

        pot_eng, enthalpy = energy_data[temperature]

        x_ref, divisor, remainder_str, n_ref, n_rem, g, comp_min_str = endmember_decompose_and_divisor(
            compound, reference_formula
        )
        if np.isnan(divisor) or divisor == 0:
            print(f"警告: {compound} divisor 无效，跳过")
            continue

        pot_norm = pot_eng / divisor - x_ref * P1 -(1 - x_ref) * P0 if not np.isnan(pot_eng) else np.nan
        ent_norm = enthalpy / divisor - x_ref * ent1 -(1 - x_ref) * ent0 if not np.isnan(enthalpy) else np.nan

        print(
            f"{compound} -> {comp_min_str} (g={g}): n_ref={n_ref}, n_rem={n_rem}, "
            f"x_ref={x_ref:.4f}, divisor={divisor}, pot_norm={pot_norm:.10f}, ent_norm={ent_norm:.10f}"
        )

        if not np.isnan(pot_norm):
            data_pot.append((x_ref, pot_norm))
        if not np.isnan(ent_norm):
            data_ent.append((x_ref, ent_norm))

    if not data_pot and not data_ent:
        print("错误: 没有有效数据可以绘图/拟合")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 势能
    if data_pot:
        data_pot.sort(key=lambda t: t[0])
        x_pot = np.array([d[0] for d in data_pot], dtype=float)
        y_pot = np.array([d[1] for d in data_pot], dtype=float)
        plot_one_panel(
            axes[0], x_pot, y_pot,
            color_points="tab:blue",
            color_line_default="red",
            title=f'Potential Energy vs x({reference_formula}) at {temperature} K',
            ylabel='Normalized Potential Energy (eV)',
            reference_formula=reference_formula
        )
    else:
        axes[0].text(0.5, 0.5, 'No valid potential energy data',
                     ha='center', va='center', transform=axes[0].transAxes)

    # 焓
    if data_ent:
        data_ent.sort(key=lambda t: t[0])
        x_ent = np.array([d[0] for d in data_ent], dtype=float)
        y_ent = np.array([d[1] for d in data_ent], dtype=float)
        plot_one_panel(
            axes[1], x_ent, y_ent,
            color_points="green",
            color_line_default="red",
            title=f'Enthalpy vs x({reference_formula}) at {temperature} K and 0 GPa',
            ylabel='Normalized Enthalpy (eV)',
            reference_formula=reference_formula
        )
    else:
        axes[1].text(0.5, 0.5, 'No valid enthalpy data',
                     ha='center', va='center', transform=axes[1].transAxes)

    plt.tight_layout()
    output_file = f"{output_prefix}_{reference_formula}_{temperature}K.eps"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n图形已保存到: {output_file}")
    plt.show()


def main():
    print(f"\n{'='*60}")
    print(f"REFERENCE_FORMULA: {REFERENCE_FORMULA}")
    print(f"TEMPERATURE: {TEMPERATURE} K")
    print(f"FILES: {len(FILES)}")
    print(f"{'='*60}\n")
    plot_energy_composition(FILES, REFERENCE_FORMULA, TEMPERATURE, OUTPUT_PREFIX)


if __name__ == '__main__':
    main()