#!/usr/bin/env python3
"""Solve an ilp_baseline .lp file with HiGHS: solve_lp_highs.py FILE.lp TIME [--relax]"""
import sys, time, highspy
h = highspy.Highs(); h.setOptionValue("output_flag", False)
h.setOptionValue("time_limit", float(sys.argv[2])); h.setOptionValue("threads", 2)
h.readModel(sys.argv[1])
if "--relax" in sys.argv[3:]:
    lp = h.getLp(); lp.integrality_ = [highspy.HighsVarType.kContinuous] * lp.num_col_; h.passModel(lp)
t = time.time(); h.run(); info = h.getInfo()
st = h.modelStatusToString(h.getModelStatus())
best = round(info.objective_function_value) if info.primal_solution_status == 2 else ""
print(f"status={st.replace(' ', '_')} best={best} bound={info.mip_dual_bound:.3f} secs={time.time()-t:.2f}")
