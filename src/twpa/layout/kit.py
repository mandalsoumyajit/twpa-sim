"""Parameterized primitives adapted from the user-supplied NY CREATES
research KIT layout generator (generate_v2_gds.py). This module preserves
its straight-comb and constant-width radial-finger construction. The biased
sheet-squares match is an approximation; it does not establish impedance.
"""
import klayout.db as db
import numpy as np
from scipy.optimize import brentq

def primitives(ly, design):
    """Straight comb, constant-width-finger U-turns and tapered launches."""
    p_straight = design.pitch
    w_line_straight = design.spine_width
    s_line = design.spine_gap
    wf = design.finger_width
    sf = design.finger_gap
    lf = design.finger_length
    llf = design.loaded_finger_length
    N_sc_cells = design.unloaded_cells + design.loaded_cells
    R_turn = design.turn_radius
    N_turn_sc = design.turn_periods
    N_turn_cells = N_turn_sc * N_sc_cells
    p_turn = np.pi * R_turn / N_turn_cells

    sc_polys = []
    y_gnd_ext = R_turn
    for c in range(N_sc_cells):
        is_loaded = (design.unloaded_cells//2 <= c < design.unloaded_cells//2 + design.loaded_cells)
        l_curr = llf if is_loaded else lf
        x0 = c * p_straight
        x1 = (c + 1) * p_straight
        xm = 0.5 * (x0 + x1)

        y_gt_bot = 0.5 * w_line_straight + s_line
        y_slot_top = 0.5 * w_line_straight + lf + s_line
        pts_gt = [
            db.DPoint(x0, y_gt_bot),
            db.DPoint(x0, y_gnd_ext),
            db.DPoint(x1, y_gnd_ext),
            db.DPoint(x1, y_gt_bot),
            db.DPoint(xm + 0.5*wf + sf, y_gt_bot),
            db.DPoint(xm + 0.5*wf + sf, y_slot_top),
            db.DPoint(xm - 0.5*wf - sf, y_slot_top),
            db.DPoint(xm - 0.5*wf - sf, y_gt_bot),
        ]
        sc_polys.append(db.DPolygon(pts_gt).to_itype(ly.dbu))

        y_gb_top = -0.5 * w_line_straight - s_line
        y_slot_bot = -0.5 * w_line_straight - lf - s_line
        pts_gb = [
            db.DPoint(x0, -y_gnd_ext),
            db.DPoint(x0, y_gb_top),
            db.DPoint(xm - 0.5*wf - sf, y_gb_top),
            db.DPoint(xm - 0.5*wf - sf, y_slot_bot),
            db.DPoint(xm + 0.5*wf + sf, y_slot_bot),
            db.DPoint(xm + 0.5*wf + sf, y_gb_top),
            db.DPoint(x1, y_gb_top),
            db.DPoint(x1, -y_gnd_ext),
        ]
        sc_polys.append(db.DPolygon(pts_gb).to_itype(ly.dbu))

        y_cl_top = 0.5 * w_line_straight
        y_cl_bot = -0.5 * w_line_straight
        y_ft_top = y_cl_top + l_curr
        y_fb_bot = y_cl_bot - l_curr
        pts_sig = [
            db.DPoint(xm - 0.5*wf, y_fb_bot),
            db.DPoint(xm - 0.5*wf, y_cl_bot),
            db.DPoint(x0, y_cl_bot),
            db.DPoint(x0, y_cl_top),
            db.DPoint(xm - 0.5*wf, y_cl_top),
            db.DPoint(xm - 0.5*wf, y_ft_top),
            db.DPoint(xm + 0.5*wf, y_ft_top),
            db.DPoint(xm + 0.5*wf, y_cl_top),
            db.DPoint(x1, y_cl_top),
            db.DPoint(x1, y_cl_bot),
            db.DPoint(xm + 0.5*wf, y_cl_bot),
            db.DPoint(xm + 0.5*wf, y_fb_bot),
        ]
        sc_polys.append(db.DPolygon(pts_sig).to_itype(ly.dbu))

    dtheta = np.pi / N_turn_cells
    th_k = np.pi/2 - (np.arange(N_turn_cells) + 0.5) * dtheta

    def w_line_func(th, w_max):
        d_end = np.minimum(np.pi/2 - th, th + np.pi/2)
        th_taper = design.taper_cells * dtheta
        alpha = np.clip(d_end / th_taper, 0.0, 1.0)
        s = alpha * alpha * (3.0 - 2.0 * alpha)
        return w_line_straight + (w_max - w_line_straight) * s

    target_sq_biased = N_turn_cells * (p_straight / w_line_straight) * (1.0 + design.bias_ratio**2)
    def turn_biased_sq_err(wm):
        w_arr = w_line_func(th_k, wm)
        sq_b = (p_turn / w_arr) * (1.0 + design.bias_ratio**2 * (w_line_straight / w_arr)**2)
        return np.sum(sq_b) - target_sq_biased

    w_max_opt = brentq(turn_biased_sq_err, w_line_straight/10, 10*w_line_straight)

    n_arc = design.arc_points
    thetas = np.linspace(-np.pi/2, np.pi/2, n_arc)
    w_th = w_line_func(thetas, w_max_opt)
    r_in_line = R_turn - 0.5 * w_th
    r_out_line = R_turn + 0.5 * w_th

    reg_sig_turn = db.Region(db.DPolygon(
        [db.DPoint(r*np.cos(t), r*np.sin(t)) for r, t in zip(r_out_line, thetas)] +
        [db.DPoint(r*np.cos(t), r*np.sin(t)) for r, t in zip(reversed(r_in_line), reversed(thetas))]
    ).to_itype(ly.dbu))

    r_gnd_in = r_in_line - s_line
    reg_gnd_in = db.Region(db.DPolygon(
        [db.DPoint(0, -r_gnd_in[0])] +
        [db.DPoint(r*np.cos(t), r*np.sin(t)) for r, t in zip(r_gnd_in, thetas)] +
        [db.DPoint(0, r_gnd_in[-1])]
    ).to_itype(ly.dbu))

    r_gnd_out_inner = r_out_line + s_line
    r_outer_box = max(design.launch_length+design.taper_length, R_turn+lf+10)
    y_outer_box = 2.0 * R_turn
    reg_gnd_out = db.Region(db.DPolygon(
        [db.DPoint(0, r_gnd_out_inner[-1])] +
        [db.DPoint(r*np.cos(t), r*np.sin(t)) for r, t in zip(reversed(r_gnd_out_inner), reversed(thetas))] +
        [db.DPoint(0, -r_gnd_out_inner[0]),
         db.DPoint(0, -y_outer_box),
         db.DPoint(r_outer_box, -y_outer_box),
         db.DPoint(r_outer_box, y_outer_box),
         db.DPoint(0, y_outer_box)]
    ).to_itype(ly.dbu))

    reg_slots_in = db.Region()
    reg_slots_out = db.Region()
    for k in range(N_turn_cells):
        is_loaded = (design.unloaded_cells//2 <= (k % N_sc_cells) < design.unloaded_cells//2 + design.loaded_cells)
        l_curr = llf if is_loaded else lf
        th = th_k[k]
        w_loc = w_line_func(th, w_max_opt)
        r_in_k = R_turn - 0.5 * w_loc
        r_out_k = R_turn + 0.5 * w_loc
        ux, uy = np.cos(th), np.sin(th)
        tx, ty = -np.sin(th), np.cos(th)

        def rect_poly(u0, u1, half_w):
            return db.DPolygon([
                db.DPoint(u0*ux - half_w*tx, u0*uy - half_w*ty),
                db.DPoint(u1*ux - half_w*tx, u1*uy - half_w*ty),
                db.DPoint(u1*ux + half_w*tx, u1*uy + half_w*ty),
                db.DPoint(u0*ux + half_w*tx, u0*uy + half_w*ty),
            ]).to_itype(ly.dbu)

        reg_sig_turn.insert(rect_poly(r_out_k - 0.1, r_out_k + l_curr, 0.5 * wf))
        reg_sig_turn.insert(rect_poly(r_in_k - l_curr, r_in_k + 0.1, 0.5 * wf))
        reg_slots_out.insert(rect_poly(R_turn, r_out_k + lf + s_line, 0.5 * wf + sf))
        reg_slots_in.insert(rect_poly(r_in_k - lf - s_line, R_turn, 0.5 * wf + sf))

    reg_turn_proto = reg_sig_turn.merged() + (reg_gnd_in - reg_slots_in).merged() + (reg_gnd_out - reg_slots_out).merged()
    turn_polys_right = list(reg_turn_proto.each())
    tr_mirror_x = db.Trans(db.Trans.M90)
    turn_polys_left = [p.transformed(tr_mirror_x) for p in turn_polys_right]

    W_pad = design.launch_width
    G_pad = design.launch_gap
    y_corr = W_pad/2 + G_pad
    L_bay = design.launch_length + design.taper_length
    L_pad = design.launch_length
    L_tap = design.taper_length

    pts_launch_sig = [
        db.DPoint(-L_bay, -0.5 * W_pad),
        db.DPoint(-L_tap, -0.5 * W_pad),
        db.DPoint(0.0, -0.5 * w_line_straight),
        db.DPoint(0.0, +0.5 * w_line_straight),
        db.DPoint(-L_tap, +0.5 * W_pad),
        db.DPoint(-L_bay, +0.5 * W_pad),
    ]
    poly_launch_sig = db.DPolygon(pts_launch_sig).to_itype(ly.dbu)

    y_gt0 = 0.5 * w_line_straight + s_line
    pts_launch_gt = [
        db.DPoint(-L_bay, y_corr),
        db.DPoint(-L_bay, R_turn),
        db.DPoint(0.0, R_turn),
        db.DPoint(0.0, y_gt0),
        db.DPoint(-design.finger_gap, y_gt0),
        db.DPoint(-L_tap, y_corr),
    ]
    poly_launch_gt = db.DPolygon(pts_launch_gt).to_itype(ly.dbu)

    pts_launch_gb = [
        db.DPoint(-L_bay, -R_turn),
        db.DPoint(-L_bay, -y_corr),
        db.DPoint(-L_tap, -y_corr),
        db.DPoint(-design.finger_gap, -y_gt0),
        db.DPoint(0.0, -y_gt0),
        db.DPoint(0.0, -R_turn),
    ]
    poly_launch_gb = db.DPolygon(pts_launch_gb).to_itype(ly.dbu)

    launch_polys_left = [poly_launch_sig, poly_launch_gt, poly_launch_gb]
    launch_polys_right = [p.transformed(tr_mirror_x) for p in launch_polys_left]

    return sc_polys, turn_polys_right, turn_polys_left, launch_polys_left, launch_polys_right, p_straight, N_sc_cells, R_turn, w_max_opt

