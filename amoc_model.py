"""
================================================================================
  Rooth 3-Box AMOC Model  —  Presentation Dashboard
================================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.widgets import Slider, Button
import matplotlib.patheffects as pe

# ──────────────────────────────────────────────────────────────────────────────
#  PHYSICS  
# ──────────────────────────────────────────────────────────────────────────────

RHO_0  = 1025.0
ALPHA  = 1.5e-4
BETA   = 8.0e-4
T_REF  = 15.0
S_REF  = 35.0
SV     = 1.0e6

V_N = 1.7e16
V_T = 8.0e16
V_S = 1.7e16
K   = 4.0e7

T_N_ATM = 5.0
T_T_ATM = 25.0
T_S_ATM = 2.0
TAU_T   = 30 * 86_400
Q_WIND_BASE = 4.0e6

T_YEARS   = 200
T_EVAL_YR = np.linspace(0, T_YEARS, 2000)
T_EVAL    = T_EVAL_YR * 365 * 86_400

Y0 = [5.0, 24.0, 3.0, 35.5, 36.5, 34.5]   # [T_N, T_T, T_S, S_N, S_T, S_S]


def density(T, S):
    return RHO_0 * (-ALPHA * (T - T_REF) + BETA * (S - S_REF))


def amoc_flow(T_N, T_S, S_N, S_S):
    return K * (density(T_N, S_N) - density(T_S, S_S))


def rooth_3box(t, y, params):
    T_N, T_T, T_S, S_N, S_T, S_S = y
    t_yr  = t / (365 * 86_400)
    ramp  = max(params['ramp_years'], 1.0)
    scale = min(t_yr / ramp, 1.0)

    F_grl  = params['F_greenland'] * scale * SV
    F_ant  = params['F_antarctic'] * scale * SV
    dT_atm = params['delta_T_atm'] * scale
    w_fac  = params['wind_factor']

    Tn_atm = T_N_ATM + dT_atm * 0.70
    Tt_atm = T_T_ATM + dT_atm * 1.00
    Ts_atm = T_S_ATM + dT_atm * 0.55

    q   = amoc_flow(T_N, T_S, S_N, S_S)
    q_w = w_fac * Q_WIND_BASE

    if q >= 0:
        adv_TN = q*(T_T-T_N)/V_N;  adv_TT = q*(T_S-T_T)/V_T;  adv_TS = q*(T_N-T_S)/V_S
        adv_SN = q*(S_T-S_N)/V_N;  adv_ST = q*(S_S-S_T)/V_T;  adv_SS = q*(S_N-S_S)/V_S
    else:
        qa = abs(q)
        adv_TN = qa*(T_S-T_N)/V_N; adv_TT = qa*(T_N-T_T)/V_T; adv_TS = qa*(T_T-T_S)/V_S
        adv_SN = qa*(S_S-S_N)/V_N; adv_ST = qa*(S_N-S_T)/V_T; adv_SS = qa*(S_T-S_S)/V_S

    wd_TN = q_w*(T_T-T_N)/V_N;  wd_TS = q_w*(T_T-T_S)/V_S
    wd_SN = q_w*(S_T-S_N)/V_N;  wd_SS = q_w*(S_T-S_S)/V_S

    fw_SN = -F_grl*S_REF/V_N
    fw_SS = -F_ant*S_REF/V_S
    fw_ST =  (F_grl+F_ant)*S_REF/V_T

    dTN = adv_TN + wd_TN + (Tn_atm-T_N)/TAU_T
    dTT = adv_TT + 0      + (Tt_atm-T_T)/TAU_T
    dTS = adv_TS + wd_TS  + (Ts_atm-T_S)/TAU_T
    dSN = adv_SN + wd_SN  + fw_SN
    dST = adv_ST + 0      + fw_ST
    dSS = adv_SS + wd_SS  + fw_SS

    return [dTN, dTT, dTS, dSN, dST, dSS]


def run_simulation(params):
    sol = solve_ivp(
        fun=lambda t, y: rooth_3box(t, y, params),
        t_span=(0, T_EVAL[-1]), y0=Y0, t_eval=T_EVAL,
        method='RK45', rtol=1e-7, atol=1e-10,
    )
    T_N, T_T, T_S = sol.y[0], sol.y[1], sol.y[2]
    S_N, S_T, S_S = sol.y[3], sol.y[4], sol.y[5]
    q_Sv = np.array([amoc_flow(T_N[i], T_S[i], S_N[i], S_S[i])/SV
                     for i in range(len(sol.t))])
    return dict(years=T_EVAL_YR,
                T_N=T_N, T_T=T_T, T_S=T_S,
                S_N=S_N, S_T=S_T, S_S=S_S, q_Sv=q_Sv)


# ──────────────────────────────────────────────────────────────────────────────
#  COLOR PALETTE
# ──────────────────────────────────────────────────────────────────────────────

BG       = '#0d1117'
PANEL    = '#161b22'
BORDER   = '#30363d'
TXT      = '#e6edf3'
TXT_DIM  = '#8b949e'
BLUE     = '#58a6ff'
ORANGE   = '#f0883e'
GREEN    = '#3fb950'
RED      = '#f85149'
YELLOW   = '#d29922'
PURPLE   = '#bc8cff'
TEAL     = '#39d353'

C_NORTH  = '#58a6ff'   # North box colour
C_TROP   = '#f0883e'   # Tropics box colour
C_SOUTH  = '#39d353'   # South box colour


# ──────────────────────────────────────────────────────────────────────────────
#  SCHEMATIC  — draw the 3-box cartoon
# ──────────────────────────────────────────────────────────────────────────────

def draw_schematic(ax):
    """Draw an annotated 3-box diagram of the AMOC loop."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.set_facecolor(PANEL)
    ax.axis('off')
    ax.set_title('How the 3-Box Model Works', color=TXT, fontsize=10, fontweight='bold', pad=6)

    # ── Boxes ────────────────────────────────────────────────────────────────
    box_style = dict(boxstyle='round,pad=0.3', lw=1.5)

    boxes = [
        # (label, sublabel, x, y, w, h, color)
        ('NORTH\nAtlantic',  'Cold & salty\n→ dense, sinks',  0.4, 2.8, 2.2, 2.4, C_NORTH),
        ('TROPICS',          'Warm & salty\n', 3.9, 2.8, 2.2, 2.4, C_TROP),
        ('SOUTH\nAtlantic',  'Cold & fresh\n→ lighter',        7.4, 2.8, 2.2, 2.4, C_SOUTH),
    ]

    for label, sub, x, y, w, h, col in boxes:
        rect = FancyBboxPatch((x, y), w, h,
                              boxstyle='round,pad=0.15',
                              facecolor=col+'22', edgecolor=col, lw=1.8,
                              transform=ax.transData)
        ax.add_patch(rect)
        ax.text(x+w/2, y+h*0.62, label, ha='center', va='center',
                color=col, fontsize=8.5, fontweight='bold')
        ax.text(x+w/2, y+h*0.25, sub, ha='center', va='center',
                color=TXT_DIM, fontsize=6.5)

    # ── Arrows: surface flow (warm, top) ────────────────────────────────────
    arrow_kw = dict(arrowstyle='->', color=ORANGE, lw=2,
                    mutation_scale=14, transform=ax.transData)

    # Surface: S → T  (warm surface flow northward, arrow points toward tropics)
    ax.annotate('', xy=(6.1, 5.5), xytext=(7.4, 5.5),
                arrowprops=dict(arrowstyle='->', color=ORANGE, lw=2, mutation_scale=14))
    ax.text(6.75, 5.7, 'warm surface water\nflows north', ha='center',
            color=ORANGE, fontsize=6)

    # Surface: T → N  (warm surface flow continues northward)
    ax.annotate('', xy=(2.6, 5.5), xytext=(3.9, 5.5),
                arrowprops=dict(arrowstyle='->', color=ORANGE, lw=2, mutation_scale=14))
    ax.text(3.25, 5.7, 'warm water\nheats the north', ha='center',
            color=ORANGE, fontsize=6)

    # ── Sinking in North ─────────────────────────────────────────────────────
    ax.annotate('', xy=(1.5, 2.8), xytext=(1.5, 5.2),
                arrowprops=dict(arrowstyle='->', color=BLUE, lw=2, mutation_scale=14))
    ax.text(0.05, 4.0, 'dense\nwater\nsinks', ha='center',
            color=BLUE, fontsize=6, rotation=90)

    # ── Deep flow: N → S ─────────────────────────────────────────────────────
    ax.annotate('', xy=(8.5, 1.5), xytext=(1.5, 1.5),
                arrowprops=dict(arrowstyle='->', color=BLUE, lw=2, mutation_scale=14))
    ax.text(5.0, 1.2, 'cold deep water flows south along ocean floor',
            ha='center', color=BLUE, fontsize=6)

    # ── Upwelling in South ───────────────────────────────────────────────────
    ax.annotate('', xy=(8.5, 5.2), xytext=(8.5, 2.8),
                arrowprops=dict(arrowstyle='->', color=C_SOUTH, lw=2, mutation_scale=14))
    ax.text(9.95, 4.0, 'water\nrises\nback up', ha='center',
            color=C_SOUTH, fontsize=6, rotation=90)

    # ── Forcing labels ───────────────────────────────────────────────────────
    ax.text(1.5, 0.55, '← Greenland melt adds\nfreshwater here',
            ha='center', color=C_NORTH, fontsize=6.5,
            bbox=dict(fc=C_NORTH+'18', ec=C_NORTH, lw=0.8, boxstyle='round,pad=0.2'))
    ax.text(8.5, 0.55, 'Antarctic melt adds\nfreshwater here →',
            ha='center', color=C_SOUTH, fontsize=6.5,
            bbox=dict(fc=C_SOUTH+'18', ec=C_SOUTH, lw=0.8, boxstyle='round,pad=0.2'))
    ax.text(5.0, 0.1, '▲ The AMOC carries ~1.3 petawatts of heat to Europe',
            ha='center', color=TXT_DIM, fontsize=6, style='italic')


# ──────────────────────────────────────────────────────────────────────────────
#  ANNOTATED GRAPH HELPERS
# ──────────────────────────────────────────────────────────────────────────────

def style_ax(ax, title, xlabel, ylabel, caption=''):
    """Apply consistent dark styling + optional caption."""
    ax.set_facecolor(PANEL)
    ax.set_title(title, color=TXT, fontsize=9.5, fontweight='bold', pad=5)
    ax.set_xlabel(xlabel, color=TXT_DIM, fontsize=8)
    ax.set_ylabel(ylabel, color=TXT_DIM, fontsize=8)
    ax.tick_params(colors=TXT_DIM, labelsize=7.5)
    for sp in ax.spines.values():
        sp.set_edgecolor(BORDER)
    ax.grid(True, color=BORDER, alpha=0.5, linewidth=0.6)
    ax.set_xlim(0, T_YEARS)
    if caption:
        ax.text(0.5, -0.28, caption, transform=ax.transAxes,
                ha='center', va='top', color=TXT_DIM, fontsize=6.8, style='italic',
                wrap=True)


def add_threshold(ax, y_val, label, color, side='right'):
    ax.axhline(y_val, color=color, lw=1, ls='--', alpha=0.7)
    x_pos = T_YEARS * 0.97 if side == 'right' else T_YEARS * 0.03
    ha = 'right' if side == 'right' else 'left'
    ax.text(x_pos, y_val, f' {label}', ha=ha, va='bottom',
            color=color, fontsize=7, alpha=0.9)


# ──────────────────────────────────────────────────────────────────────────────
#  MAIN DASHBOARD
# ──────────────────────────────────────────────────────────────────────────────

def launch_dashboard():
    DEFAULT = dict(F_greenland=0.0, F_antarctic=0.0,
                   delta_T_atm=0.0, wind_factor=1.0, ramp_years=30)

    res = run_simulation(DEFAULT)

    # ── Figure ───────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(17, 10))
    fig.patch.set_facecolor(BG)
    fig.suptitle(
        'AMOC 3-Box Model',
        color=TXT, fontsize=14, fontweight='bold', y=0.98
    )

    # Grid: 3 rows × 4 cols, right column for sliders
    gs = gridspec.GridSpec(
        3, 4,
        figure=fig,
        left=0.05, right=0.70,
        top=0.93, bottom=0.10,
        hspace=0.60, wspace=0.38,
    )

    ax_schema = fig.add_subplot(gs[0, :2])   # schematic — top left
    ax_q      = fig.add_subplot(gs[0, 2:])   # AMOC strength — top right
    ax_sal    = fig.add_subplot(gs[1, :2])   # North salinity
    ax_temp   = fig.add_subplot(gs[1, 2:])   # Box temperatures
    ax_status = fig.add_subplot(gs[2, :2])   # Status / explanation
    ax_rho    = fig.add_subplot(gs[2, 2:])   # Density driver

    # ── Schematic ────────────────────────────────────────────────────────────
    draw_schematic(ax_schema)

    # ── AMOC strength ────────────────────────────────────────────────────────
    ln_q, = ax_q.plot(res['years'], res['q_Sv'], color=BLUE, lw=2.2, zorder=5)
    ax_q.fill_between(res['years'], res['q_Sv'], 0,
                      where=res['q_Sv'] > 0, alpha=0.12, color=BLUE, zorder=3)
    ax_q.fill_between(res['years'], res['q_Sv'], 0,
                      where=res['q_Sv'] < 0, alpha=0.18, color=RED, zorder=3)

    add_threshold(ax_q, 17,  'Today ~17 Sv', GREEN)
    add_threshold(ax_q,  5,  'Severely weakened', YELLOW)
    add_threshold(ax_q,  0,  'Collapse / reversal', RED)

    style_ax(ax_q,
             'AMOC Overturning Strength  (the main output)',
             'Year', 'Flow rate  (Sv = million m³/s)',
             'Higher = stronger circulation, more heat reaches Europe.\n'
             '0 = AMOC has stopped.  Negative = circulation has reversed.')

    # ── North Atlantic Salinity ───────────────────────────────────────────────
    ln_SN, = ax_sal.plot(res['years'], res['S_N'], color=C_NORTH, lw=2, label='North Atlantic')
    ln_ST, = ax_sal.plot(res['years'], res['S_T'], color=C_TROP,  lw=2, label='Tropics', ls='--')
    ln_SS, = ax_sal.plot(res['years'], res['S_S'], color=C_SOUTH, lw=2, label='South Atlantic', ls=':')

    ax_sal.legend(fontsize=7.5, facecolor=PANEL, labelcolor=TXT,
                  framealpha=0.8, loc='upper right')
    style_ax(ax_sal,
             'Salinity in Each Box  (salt concentration)',
             'Year', 'Salinity  (psu — practical salinity units)',
             'Freshwater from melting ice dilutes the North Atlantic (blue drops).\n'
             'Less salt = less dense = water stops sinking = AMOC weakens.')

    add_threshold(ax_sal, 35.0, '35 psu reference', TXT_DIM, side='left')

    # ── Box temperatures ──────────────────────────────────────────────────────
    ln_TN, = ax_temp.plot(res['years'], res['T_N'], color=C_NORTH, lw=2, label='North Atlantic')
    ln_TT, = ax_temp.plot(res['years'], res['T_T'], color=C_TROP,  lw=2, label='Tropics', ls='--')
    ln_TS, = ax_temp.plot(res['years'], res['T_S'], color=C_SOUTH, lw=2, label='South Atlantic', ls=':')

    ax_temp.legend(fontsize=7.5, facecolor=PANEL, labelcolor=TXT,
                   framealpha=0.8, loc='upper right')
    style_ax(ax_temp,
             'Temperature in Each Box  (ocean heat)',
             'Year', 'Temperature  (°C)',
             'If AMOC weakens, less warm water reaches the North — it cools.\n'
             'The tropics may warm further as heat is no longer exported north.')

    # ── Density driver ────────────────────────────────────────────────────────
    rho_N = density(res['T_N'], res['S_N'])
    rho_S = density(res['T_S'], res['S_S'])
    ln_dr, = ax_rho.plot(res['years'], rho_N - rho_S, color=PURPLE, lw=2)
    ax_rho.fill_between(res['years'], rho_N - rho_S, 0,
                        where=(rho_N - rho_S) > 0, alpha=0.12, color=PURPLE)
    ax_rho.fill_between(res['years'], rho_N - rho_S, 0,
                        where=(rho_N - rho_S) < 0, alpha=0.18, color=RED)
    add_threshold(ax_rho, 0, 'AMOC stops here', RED, side='left')

    style_ax(ax_rho,
             'Density Difference  ρ_North − ρ_South  (the engine)',
             'Year', 'Δρ  (kg/m³)',
             'This is literally what drives the AMOC.\n'
             'When northern water is denser than southern water, it sinks and the circulation runs.\n'
             'When this difference hits 0, the engine switches off.')

    # ── Status panel ──────────────────────────────────────────────────────────
    ax_status.set_facecolor(PANEL)
    ax_status.axis('off')
    for sp in ax_status.spines.values():
        sp.set_edgecolor(BORDER)

    ax_status.set_title('What Is Happening Right Now?', color=TXT,
                         fontsize=9.5, fontweight='bold', pad=5)

    status_main = ax_status.text(0.5, 0.72, '', transform=ax_status.transAxes,
                                  ha='center', va='center', fontsize=13,
                                  fontweight='bold', color=GREEN)

    status_detail = ax_status.text(0.5, 0.38, '', transform=ax_status.transAxes,
                                    ha='center', va='center', fontsize=8,
                                    color=TXT_DIM, wrap=True,
                                    multialignment='center')

    status_vals = ax_status.text(0.5, 0.08, '', transform=ax_status.transAxes,
                                  ha='center', va='center', fontsize=7.5,
                                  color=TXT_DIM, family='monospace')

    # ── Sliders ───────────────────────────────────────────────────────────────
    SL_LEFT  = 0.725
    SL_WIDTH = 0.240

    slider_defs = [
        # (y_bottom, label, min, max, init, step, unit_note)
        (0.820, 'Greenland Ice Melt  (freshwater into North Atlantic)',
         0.0, 0.50, 0.0, 0.005,
         '0 = no extra melt     0.5 Sv = major collapse scenario'),

        (0.675, 'Antarctic Ice Melt  (freshwater into South Atlantic)',
         0.0, 0.30, 0.0, 0.005,
         '0 = no extra melt     0.3 Sv = significant melt'),

        (0.530, 'Atmospheric Warming  (above pre-industrial levels)',
         0.0, 6.0, 0.0, 0.1,
         '0 = pre-industrial     2°C = Paris target     6°C = extreme scenario'),

        (0.385, 'Trade Wind Strength  (wind-driven ocean mixing)',
         0.0, 2.0, 1.0, 0.05,
         '0 = no wind     1 = present-day     2 = doubled trade winds'),

        (0.240, 'How quickly forcing builds up  (ramp-up years)',
         5.0, 100.0, 30.0, 5.0,
         '5 = abrupt change     50 = gradual change over decades'),
    ]

    sliders = []
    for y_bot, label, lo, hi, init, step, note in slider_defs:
        # Label above slider
        fig.text(SL_LEFT, y_bot + 0.062, label,
                 color=TXT, fontsize=8, fontweight='bold')
        fig.text(SL_LEFT, y_bot + 0.043, note,
                 color=TXT_DIM, fontsize=6.5, style='italic')

        ax_sl = fig.add_axes([SL_LEFT, y_bot, SL_WIDTH, 0.025])
        ax_sl.set_facecolor(PANEL)
        sl = Slider(ax=ax_sl, label='', valmin=lo, valmax=hi,
                    valinit=init, valstep=step,
                    color=PURPLE, track_color=BORDER)
        sl.valtext.set_color(PURPLE)
        sl.valtext.set_fontsize(8)
        sliders.append(sl)

    sl_grl, sl_ant, sl_dT, sl_wnd, sl_rmp = sliders

    # Section header for sliders
    fig.text(SL_LEFT + SL_WIDTH/2, 0.945,
             'Adjust Forcing Parameters',
             ha='center', color=TXT, fontsize=11, fontweight='bold')
    fig.text(SL_LEFT + SL_WIDTH/2, 0.925,
             'Drag any slider to see how each factor affects the AMOC over 200 years',
             ha='center', color=TXT_DIM, fontsize=7.5, style='italic')

    # Reset button
    ax_btn = fig.add_axes([SL_LEFT + SL_WIDTH/2 - 0.05, 0.115, 0.10, 0.035])
    btn = Button(ax_btn, '↺  Reset All', color='#21262d', hovercolor=PURPLE)
    btn.label.set_color(TXT)
    btn.label.set_fontsize(8.5)

    # ── Update ────────────────────────────────────────────────────────────────
    def update(_=None):
        p = dict(
            F_greenland = sl_grl.val,
            F_antarctic = sl_ant.val,
            delta_T_atm = sl_dT.val,
            wind_factor = sl_wnd.val,
            ramp_years  = sl_rmp.val,
        )
        r = run_simulation(p)
        yr = r['years']

        # ── Update lines ──────────────────────────────────────────────────
        ln_q.set_ydata(r['q_Sv'])
        for c in ax_q.collections[:]: c.remove()
        ax_q.fill_between(yr, r['q_Sv'], 0, where=r['q_Sv'] > 0,
                           alpha=0.12, color=BLUE, zorder=3)
        ax_q.fill_between(yr, r['q_Sv'], 0, where=r['q_Sv'] < 0,
                           alpha=0.18, color=RED, zorder=3)
        ax_q.relim(); ax_q.autoscale_view()

        ln_SN.set_ydata(r['S_N']); ln_ST.set_ydata(r['S_T']); ln_SS.set_ydata(r['S_S'])
        ax_sal.relim(); ax_sal.autoscale_view()

        ln_TN.set_ydata(r['T_N']); ln_TT.set_ydata(r['T_T']); ln_TS.set_ydata(r['T_S'])
        ax_temp.relim(); ax_temp.autoscale_view()

        rn = density(r['T_N'], r['S_N'])
        rs = density(r['T_S'], r['S_S'])
        dr = rn - rs
        ln_dr.set_ydata(dr)
        for c in ax_rho.collections[:]: c.remove()
        ax_rho.fill_between(yr, dr, 0, where=dr > 0, alpha=0.12, color=PURPLE)
        ax_rho.fill_between(yr, dr, 0, where=dr < 0, alpha=0.18, color=RED)
        ax_rho.relim(); ax_rho.autoscale_view()

        # ── Status text ───────────────────────────────────────────────────
        q_now   = r['q_Sv'][-1]
        q_min   = r['q_Sv'].min()
        q_start = r['q_Sv'][0]
        q_change = q_now - q_start
        dS = r['S_N'][-1] - r['S_N'][0]
        dT = r['T_N'][-1] - r['T_N'][0]

        if q_min < 0:
            headline = '⚠  AMOC HAS COLLAPSED & REVERSED'
            detail = (
                f'The circulation has completely shut down and started running backwards.\n'
                f'This would mean dramatically colder winters across northern Europe,\n'
                f'faster sea-level rise on the US East Coast, and disrupted monsoons globally.\n'
                f'Final strength: {q_now:.1f} Sv   (minimum reached: {q_min:.1f} Sv)'
            )
            col = RED
        elif q_now < 5:
            headline = f'⚡  AMOC SEVERELY WEAKENED  ({q_now:.1f} Sv)'
            detail = (
                f'The AMOC is operating at critically low levels — down {abs(q_change):.1f} Sv from the start.\n'
                f'At this strength, significantly less heat reaches northern Europe.\n'
                f'North Atlantic temperature has changed by {dT:+.2f}°C; '
                f'salinity changed by {dS:+.3f} psu.\n'
                f'This is near or past the point of no return.'
            )
            col = ORANGE
        elif q_now < 12:
            headline = f'⚠  AMOC NOTICEABLY WEAKENED  ({q_now:.1f} Sv)'
            detail = (
                f'The AMOC has weakened by {abs(q_change):.1f} Sv from baseline — a significant reduction.\n'
                f'Heat transport to Europe is declining. North Atlantic cooling by {dT:+.2f}°C,\n'
                f'salinity dropping by {abs(dS):.3f} psu as freshwater dilutes the northern box.\n'
                f'Continued forcing could push this past a tipping point.'
            )
            col = YELLOW
        else:
            headline = f'✓  AMOC IS HEALTHY  ({q_now:.1f} Sv)'
            detail = (
                f'The circulation is running at {q_now:.1f} Sv — close to the modern observed value of ~17 Sv.\n'
                f'Density difference Δρ = {dr[-1]:.3f} kg/m³ — the North is denser than the South, so water sinks.\n'
                f'North Atlantic temperature change: {dT:+.2f}°C   |   Salinity change: {dS:+.3f} psu'
            )
            col = GREEN

        status_main.set_text(headline)
        status_main.set_color(col)
        status_detail.set_text(detail)

        vals = (
            f'q_end = {q_now:.2f} Sv   '
            f'ΔS_N = {dS:+.3f} psu   '
            f'ΔT_N = {dT:+.2f}°C   '
            f'Δρ = {dr[-1]:+.4f} kg/m³'
        )
        status_vals.set_text(vals)

        fig.canvas.draw_idle()

    def reset(_):
        for sl in sliders:
            sl.reset()

    for sl in sliders:
        sl.on_changed(update)
    btn.on_clicked(reset)

    update()  # populate on first draw

    plt.show()


# ──────────────────────────────────────────────────────────────────────────────
#  ENTRY POINT
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print('=' * 68)
    print('  Rooth 3-Box AMOC')
    print('  Drag sliders to explore. Close window to exit.')
    print('=' * 68)
    launch_dashboard()