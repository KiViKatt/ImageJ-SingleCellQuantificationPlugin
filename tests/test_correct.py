# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach
# Copyright (C) 2026 Maksym Taran

# Cell Quant Studio Universal - tests of the error-corrected counts
# (cellquantcorrect.py)
#
# Loaded by run_tests.py, which provides the helpers used here.
# Every expected number is derived by hand in the comments; none is taken
# from the code under test.
# The source is kept plain ASCII with no encoding declaration on purpose.

CORRECT_TEST_EACUTE = b"\xc3\xa9".decode("utf-8")

# The worked example used by several tests: five validated fields.
#   detected D = 10, 20, 30, 40, 50     manual M = 12, 25, 31, 46, 56
CORRECT_TEST_D = [10, 20, 30, 40, 50]
CORRECT_TEST_M = [12, 25, 31, 46, 56]
# Two-sided 95 percent points of Student's t from any printed table.
CORRECT_TEST_T3 = 3.182446305
CORRECT_TEST_T4 = 2.776445105
CORRECT_TEST_T5 = 2.570581836


def correct_test_raises(function):
    try:
        function()
    except ValueError:
        return True
    return False


def correct_test_write_bytes(path, data):
    handle = open(path, "wb")
    try:
        handle.write(data)
    finally:
        handle.close()


def correct_test_read_bytes(path):
    handle = open(path, "rb")
    try:
        return handle.read()
    finally:
        handle.close()


def correct_test_column(table, name):
    index = table["headers"].index(name)
    return [row[index] for row in table["rows"]]


# -----------------------------------------------------------------------------
# 1. Exact data and the t quantile
# -----------------------------------------------------------------------------

def test_correct_exact_fits(E, t, tmp):
    """Pairs that lie exactly on a known line give that line back, with no
    scatter, so every interval collapses onto the corrected value."""
    # M = 2 + 1.1 D exactly.
    detected = [10, 20, 30, 40, 50]
    manual = [13, 24, 35, 46, 57]
    fit = E.fit_count_calibration(detected, manual, E.CORRECT_MODEL_LINE)
    t.close(fit["slope"], 1.1, "exact line: slope", 1e-9)
    t.close(fit["intercept"], 2.0, "exact line: intercept", 1e-9)
    t.check(fit["residual_sd"] < 1e-9, "exact line: no residual scatter")
    t.close(fit["r_squared"], 1.0, "exact line: R squared is 1", 1e-9)
    t.equal(fit["df"], 3, "line: n - 2 degrees of freedom")
    t.check(fit["slope_high"] - fit["slope_low"] < 1e-9, "exact line: the slope interval collapses")
    # 2 + 1.1 x 25 = 29.5
    result = E.predict_corrected_count(fit, 25)
    t.close(result["estimate"], 29.5, "exact line: corrected count of 25 detected", 1e-9)
    t.check(result["upper"] - result["lower"] < 1e-6, "exact line: the prediction interval collapses")
    t.equal((result["extrapolated"], result["clipped"]), (False, False), "inside the calibrated range, nothing flagged")

    # M = 1.25 D exactly (all values exact in binary, so the scatter is exactly 0).
    detected = [8, 16, 24, 40]
    manual = [10, 20, 30, 50]
    fit = E.fit_count_calibration(detected, manual, E.CORRECT_MODEL_ORIGIN)
    t.equal(fit["slope"], 1.25, "exact proportion: slope")
    t.equal(fit["intercept"], 0.0, "proportional model: the intercept is fixed at zero")
    t.equal(fit["residual_sd"], 0.0, "exact proportion: residual SD is exactly zero")
    t.equal(fit["df"], 3, "proportional model: n - 1 degrees of freedom")
    result = E.predict_corrected_count(fit, 32)
    t.equal((result["estimate"], result["lower"], result["upper"]), (40.0, 40.0, 40.0), "exact proportion: 32 detected is 40, no interval")
    # Through-origin slope of the line data above: sum(D M) / sum(D D)
    #   = (130 + 480 + 1050 + 1840 + 2850) / 5500 = 6350 / 5500.
    fit = E.fit_count_calibration([10, 20, 30, 40, 50], [13, 24, 35, 46, 57], E.CORRECT_MODEL_ORIGIN)
    t.close(fit["slope"], 6350.0 / 5500.0, "through-origin slope is sum(D M) / sum(D D)")

    # Student t, 97.5 percent point, against a printed table.
    printed = {1: 12.7062, 2: 4.30265, 3: 3.18245, 4: 2.77645, 5: 2.57058, 10: 2.22814,
               20: 2.08596, 30: 2.04227, 40: 2.02108, 60: 2.00030, 120: 1.97993}
    for df in sorted(printed.keys()):
        t.close(E.correct_t975(df)[0], printed[df], "t quantile for " + str(df) + " degrees of freedom", 2e-5)
        t.close(E.correct_t975_builtin(df), printed[df], "built-in t quantile for " + str(df) + " degrees of freedom", 2e-5)
    if E.CorrectTDistribution is not None:
        worst = 0.0
        for df in list(range(1, 61)) + [80, 100, 250, 1000, 100000]:
            worst = max(worst, abs(E.correct_t975(df)[0] - E.correct_t975_builtin(df)))
        t.check(worst < 1e-6, "the built-in t quantile agrees with Commons Math to 0.000001 (worst " + repr(worst) + ")")
        t.check("Commons" in E.correct_t975(5)[1], "the source of the t quantile is recorded")
    t.check(E.correct_t975(0)[0] != E.correct_t975(0)[0], "no t quantile for zero degrees of freedom")


# -----------------------------------------------------------------------------
# 2. A noisy example worked by hand
# -----------------------------------------------------------------------------

def test_correct_hand_worked_example(E, t, tmp):
    """Five noisy pairs; slope, intercept, standard errors, intervals,
    agreement and leave-one-out errors are derived by hand below."""
    D = CORRECT_TEST_D
    M = CORRECT_TEST_M

    # --- line with intercept, ordinary least squares ---
    # mean D = 30, mean M = 34
    # Sxx = 400 + 100 + 0 + 100 + 400 = 1000
    # Sxy = (-20)(-22) + (-10)(-9) + 0 + (10)(12) + (20)(22) = 1090
    # slope = 1.09, intercept = 34 - 1.09 x 30 = 1.3
    # fitted 12.2 23.1 34.0 44.9 55.8; residuals -0.2 1.9 -3.0 1.1 0.2
    # SSE = 0.04 + 3.61 + 9 + 1.21 + 0.04 = 13.9 on 3 degrees of freedom
    # residual SD = sqrt(13.9 / 3) = 2.152518
    # SE(slope) = SD / sqrt(1000) = 0.068069
    # SE(intercept) = SD x sqrt(1/5 + 900/1000) = SD x sqrt(1.1) = 2.257580
    # Syy = 484 + 81 + 9 + 144 + 484 = 1202; R squared = 1 - 13.9/1202 = 0.988436
    sd = math.sqrt(13.9 / 3.0)
    line = E.fit_count_calibration(D, M, E.CORRECT_MODEL_LINE)
    t.close(line["slope"], 1.09, "line: slope")
    t.close(line["intercept"], 1.3, "line: intercept")
    t.close(line["residual_sd"], sd, "line: residual SD")
    t.close(line["residual_sd"], 2.152518, "line: residual SD as a number", 1e-6)
    t.close(line["slope_se"], sd / math.sqrt(1000.0), "line: standard error of the slope")
    t.close(line["intercept_se"], sd * math.sqrt(1.1), "line: standard error of the intercept")
    t.close(line["t"], CORRECT_TEST_T3, "line: t quantile on 3 degrees of freedom")
    t.close(line["slope_low"], 1.09 - CORRECT_TEST_T3 * sd / math.sqrt(1000.0), "line: slope lower limit")
    t.close(line["slope_high"], 1.09 + CORRECT_TEST_T3 * sd / math.sqrt(1000.0), "line: slope upper limit")
    # intercept CI = 1.3 +/- 3.182446 x 2.257580 = -5.8846 to 8.4846 (holds zero)
    t.close(line["intercept_low"], -5.884627, "line: intercept lower limit", 1e-6)
    t.close(line["intercept_high"], 8.484627, "line: intercept upper limit", 1e-6)
    t.close(line["r_squared"], 1.0 - 13.9 / 1202.0, "line: R squared")
    # AICc = n ln(SSE/n) + 2k + 2k(k+1)/(n-k-1) with k = 3 (two coefficients
    # and the variance): 5 ln(2.78) + 6 + 24 = 35.112255
    t.close(line["aicc"], 5.0 * math.log(13.9 / 5.0) + 30.0, "line: AICc")
    # Prediction for 35 detected: 1.3 + 1.09 x 35 = 39.45
    # half width = t x SD x sqrt(1 + 1/5 + (35-30)^2/1000) = t x SD x sqrt(1.225) = 7.581862
    result = E.predict_corrected_count(line, 35)
    half = CORRECT_TEST_T3 * sd * math.sqrt(1.225)
    t.close(result["estimate"], 39.45, "line: corrected count of 35 detected")
    t.close(result["upper"] - result["estimate"], 7.581862, "line: half width of the prediction interval", 1e-6)
    t.close(result["lower"], 39.45 - half, "line: lower prediction limit")
    t.close(result["upper"], 39.45 + half, "line: upper prediction limit")
    # The interval must be the PREDICTION interval. The confidence interval of
    # the mean would have the factor sqrt(0.225) instead of sqrt(1.225).
    t.check(result["upper"] - result["estimate"] > 2.0 * CORRECT_TEST_T3 * sd * math.sqrt(0.225),
            "the interval is a prediction interval, not the narrower interval of the mean")

    # --- proportional model (through the origin) ---
    # sum D D = 5500, sum D M = 120 + 500 + 930 + 1840 + 2800 = 6190, sum M M = 6982
    # slope = 6190 / 5500 = 1.125455
    # SSE = 6982 - 6190^2 / 5500 = 15.436364 on 4 degrees of freedom; SD = 1.964457
    # SE(slope) = SD / sqrt(5500) = 0.026489
    sse = 6982.0 - 6190.0 * 6190.0 / 5500.0
    sd0 = math.sqrt(sse / 4.0)
    origin = E.fit_count_calibration(D, M, E.CORRECT_MODEL_ORIGIN)
    t.close(origin["slope"], 6190.0 / 5500.0, "proportional: slope")
    t.close(origin["residual_sd"], sd0, "proportional: residual SD")
    t.close(origin["residual_sd"], 1.964457, "proportional: residual SD as a number", 1e-6)
    t.close(origin["slope_se"], sd0 / math.sqrt(5500.0), "proportional: standard error of the slope")
    t.close(origin["t"], CORRECT_TEST_T4, "proportional: t quantile on 4 degrees of freedom")
    t.check(origin["intercept_se"] != origin["intercept_se"], "proportional: no standard error for a fixed intercept")
    t.close(origin["r_squared"], 1.0 - sse / 1202.0, "proportional: R squared about the mean")
    t.close(origin["r_squared_uncentered"], 1.0 - sse / 6982.0, "proportional: uncentered R squared")
    # AICc with k = 2: 5 ln(SSE/5) + 4 + 12/2 = 15.636440
    t.close(origin["aicc"], 5.0 * math.log(sse / 5.0) + 10.0, "proportional: AICc")
    # Prediction for 35 detected: 1.125455 x 35 = 39.390909
    # half width = t x SD x sqrt(1 + 35^2/5500) = 6.031098
    result = E.predict_corrected_count(origin, 35)
    t.close(result["estimate"], 35.0 * 6190.0 / 5500.0, "proportional: corrected count of 35 detected")
    t.close(result["upper"] - result["estimate"], CORRECT_TEST_T4 * sd0 * math.sqrt(1.0 + 1225.0 / 5500.0),
            "proportional: half width of the prediction interval")
    t.close(result["upper"] - result["estimate"], 6.031098, "proportional: half width as a number", 1e-6)

    # --- weights 1 / detected count (scatter grows with the count) ---
    # The weighted through-origin slope is the ratio of the totals:
    #   sum M / sum D = 170 / 150 = 1.133333
    # weighted SSE = sum (M - b D)^2 / D
    #   = 0.444444/10 + 5.444444/20 + 9/30 + 0.444444/40 + 0.444444/50 = 0.636667
    # SD = sqrt(0.636667 / 4) = 0.398957; SE(slope) = SD / sqrt(sum D) = 0.032575
    # variance of a new field with 30 detected = SD^2 x (30 + 30^2/150) = SD^2 x 36
    weighted = E.fit_count_calibration(D, M, E.CORRECT_MODEL_ORIGIN, E.CORRECT_WEIGHT_COUNT)
    t.close(weighted["slope"], 170.0 / 150.0, "weighted proportional: the slope is the ratio of the totals")
    t.close(weighted["residual_sd"], 0.398957, "weighted proportional: residual SD", 1e-6)
    t.close(weighted["slope_se"], 0.032575, "weighted proportional: standard error of the slope", 1e-5)
    result = E.predict_corrected_count(weighted, 30)
    t.close(result["estimate"], 34.0, "weighted proportional: corrected count of 30 detected")
    t.close(result["upper"] - result["estimate"], CORRECT_TEST_T4 * math.sqrt(0.636667 / 4.0 * 36.0),
            "weighted proportional: half width of the prediction interval", 1e-5)
    t.check(correct_test_raises(lambda: E.fit_count_calibration([0, 10, 20], [0, 11, 22], E.CORRECT_MODEL_ORIGIN, E.CORRECT_WEIGHT_COUNT)),
            "weights 1 / detected refuse a field with zero detected objects")

    # --- agreement of the raw counts ---
    # D - M = -2 -5 -1 -6 -6: mean -4; SD = sqrt((4 + 1 + 9 + 4 + 4) / 4) = sqrt(5.5)
    # limits of agreement = -4 +/- 1.96 sqrt(5.5) = -8.596607 and 0.596607
    # MAPE = 100 x (2/12 + 5/25 + 1/31 + 6/46 + 6/56) / 5 = 12.730047
    # RMSE = sqrt((4 + 25 + 1 + 36 + 36) / 5) = sqrt(20.4); MAE = 20 / 5 = 4
    agreement = E.count_agreement(D, M)
    t.close(agreement["mean_difference"], -4.0, "mean difference detected - manual")
    t.close(agreement["sd_difference"], math.sqrt(5.5), "SD of the differences")
    t.close(agreement["loa_low"], -8.596607, "lower limit of agreement", 1e-6)
    t.close(agreement["loa_high"], 0.596607, "upper limit of agreement", 1e-6)
    t.close(agreement["mape"], 100.0 * (2.0 / 12 + 5.0 / 25 + 1.0 / 31 + 6.0 / 46 + 6.0 / 56) / 5.0, "mean absolute percentage error")
    t.close(agreement["rmse"], math.sqrt(20.4), "RMSE of the raw counts")
    t.close(agreement["mae"], 4.0, "MAE of the raw counts")
    # A field with a manual count of zero cannot enter a percentage error.
    t.equal(E.count_agreement([1, 2, 3], [0, 2, 4])["mape_fields"], 2, "MAPE leaves out fields with zero manual count")

    # --- leave-one-out, proportional model ---
    # Without field i the slope is (6190 - Di Mi) / (5500 - Di Di):
    #   6070/5400, 5690/5100, 5260/4600, 4350/3900, 3390/3000
    # corrected: 11.240741 22.313725 34.304348 44.615385 56.5
    expected = [6070.0 / 5400 * 10, 5690.0 / 5100 * 20, 5260.0 / 4600 * 30, 4350.0 / 3900 * 40, 3390.0 / 3000 * 50]
    loo = E.leave_one_out_correction(D, M, E.CORRECT_MODEL_ORIGIN, E.CORRECT_WEIGHT_NONE)
    for index in range(5):
        t.close(loo["predictions"][index], expected[index], "leave-one-out corrected count of field " + str(index + 1))
    rmse = math.sqrt(sum([(expected[i] - M[i]) ** 2 for i in range(5)]) / 5.0)
    t.close(loo["rmse"], rmse, "leave-one-out RMSE")
    t.close(loo["rmse"], 2.043449, "leave-one-out RMSE as a number", 1e-6)
    # 2.043 is smaller than the raw RMSE 4.517, so correcting helps here.
    calibration = E.build_count_calibration(E.correct_pairs_from_counts(D, M))
    t.equal(calibration["model"], E.CORRECT_MODEL_ORIGIN, "the proportional model is the default")
    t.equal(calibration["weighting"], E.CORRECT_WEIGHT_NONE, "ordinary least squares is the default")
    t.equal(calibration["recommended"][E.CORRECT_MODEL_ORIGIN], True, "correction is recommended when it lowers the leave-one-out error")
    t.check(calibration["recommendation_texts"][E.CORRECT_MODEL_ORIGIN].startswith("Correction recommended"), "the recommendation is worded")
    t.check("error" not in calibration["fits"][E.CORRECT_MODEL_LINE], "the other model is fitted as well, for comparison")


# -----------------------------------------------------------------------------
# 3. Refusals, flags and warnings
# -----------------------------------------------------------------------------

def test_correct_refusals_and_flags(E, t, tmp):
    origin = E.CORRECT_MODEL_ORIGIN
    line = E.CORRECT_MODEL_LINE
    pairs = E.correct_pairs_from_counts

    # Fewer than three fields: refused, whatever the model.
    t.check(correct_test_raises(lambda: E.fit_count_calibration([10, 20], [11, 22], origin)), "two fields are refused (proportional)")
    t.check(correct_test_raises(lambda: E.fit_count_calibration([10, 20], [11, 22], line)), "two fields are refused (line)")
    t.check(correct_test_raises(lambda: E.build_count_calibration(pairs([10, 20], [11, 22]))), "two fields are refused when building a calibration")
    t.check(correct_test_raises(lambda: E.build_count_calibration([])), "no field at all is refused")
    t.check(not correct_test_raises(lambda: E.fit_count_calibration([10, 20, 30], [11, 22, 33], origin)), "three fields are accepted")
    t.check(correct_test_raises(lambda: E.fit_count_calibration([0, 0, 0], [1, 2, 3], origin)), "all-zero detected counts are refused")
    t.check(correct_test_raises(lambda: E.fit_count_calibration([5, 5, 5], [4, 5, 6], line)), "a line through fields with one detected count is refused")
    t.check(correct_test_raises(lambda: E.fit_count_calibration([10, -1, 30], [11, 22, 33], origin)), "a negative count is refused")
    t.check(correct_test_raises(lambda: E.fit_count_calibration([10, 20, 30], [11, 22], origin)), "unpaired lists are refused")
    t.check(correct_test_raises(lambda: E.build_count_calibration(pairs([10, 20, 30], [11, 22, 33]), "Cubic")), "an unknown model is refused")
    t.check(correct_test_raises(lambda: E.build_count_calibration(pairs([10, 20, 30], [11, 22, 33]), origin, "Magic")), "an unknown weighting is refused")
    # The comparison model that cannot be fitted is reported, not fatal.
    calibration = E.build_count_calibration(pairs([5, 5, 5], [4, 5, 6]), origin)
    t.check("error" in calibration["fits"][line], "a comparison model that cannot be fitted is reported as such")
    t.equal(E.correct_fit_rows(calibration, "run")[1]["Fit_Status"][:10], "NOT_FITTED", "and is marked NOT_FITTED in the fit table")

    # Extrapolation and clipping. D = 10 20 30 40, M = 5 20 35 50 lies exactly
    # on M = -10 + 1.5 D (mean D 25, mean M 27.5, Sxx 500, Sxy 750).
    fit = E.fit_count_calibration([10, 20, 30, 40], [5, 20, 35, 50], line)
    t.close(fit["slope"], 1.5, "clipping example: slope", 1e-9)
    t.close(fit["intercept"], -10.0, "clipping example: intercept", 1e-9)
    low = E.predict_corrected_count(fit, 4)       # -10 + 6 = -4, below zero
    t.equal(low["estimate"], 0.0, "a corrected count below zero is set to zero")
    t.equal(low["clipped"], True, "and flagged as clipped")
    t.equal(low["extrapolated"], True, "4 detected is below the calibrated range 10 to 40")
    t.equal((low["lower"], low["upper"]), (0.0, 0.0), "the limits of a clipped count are not negative")
    inside = E.predict_corrected_count(fit, 25)   # -10 + 37.5 = 27.5
    t.close(inside["estimate"], 27.5, "inside the range: corrected count", 1e-9)
    t.equal((inside["extrapolated"], inside["clipped"]), (False, False), "inside the range: no flag")
    t.equal(E.predict_corrected_count(fit, 10)["extrapolated"], False, "the smallest calibrated count is inside the range")
    t.equal(E.predict_corrected_count(fit, 40)["extrapolated"], False, "the largest calibrated count is inside the range")
    t.equal(E.predict_corrected_count(fit, 41)["extrapolated"], True, "41 detected is above the calibrated range")
    t.equal(E.predict_corrected_count(fit, None), None, "no count, no correction")
    t.equal(E.predict_corrected_count(fit, -1), None, "a negative count is not corrected")
    t.equal(E.predict_corrected_count(fit, float("nan")), None, "a missing count is not corrected")
    # A lower LIMIT below zero is set to zero without the clipped flag:
    # worked example, 1 detected: 1.125 +/- 2.776 x 1.964 x sqrt(1 + 1/5500).
    edge = E.predict_corrected_count(E.fit_count_calibration(CORRECT_TEST_D, CORRECT_TEST_M, origin), 1)
    t.equal((edge["lower"], edge["clipped"]), (0.0, False), "a lower limit below zero is set to zero; the count itself is not flagged")
    t.check(edge["upper"] > edge["estimate"] > 0, "the upper limit is kept")

    # Correction that does not help. D = 10..50, M = 12 18 32 38 52:
    # raw errors -2 +2 -2 +2 -2, raw RMSE = 2.
    # sum D M = 120 + 360 + 960 + 1520 + 2600 = 5560; leave-one-out slopes
    # 5440/5400 5200/5100 4600/4600 4040/3900 2960/3000 give corrected counts
    # 10.074 20.392 30 41.436 49.333, errors -1.926 2.392 -2 3.436 -2.667,
    # RMSE = 2.543545, which is LARGER than 2.
    D = CORRECT_TEST_D
    M = [12, 18, 32, 38, 52]
    expected = [5440.0 / 5400 * 10, 5200.0 / 5100 * 20, 30.0, 4040.0 / 3900 * 40, 2960.0 / 3000 * 50]
    rmse = math.sqrt(sum([(expected[i] - M[i]) ** 2 for i in range(5)]) / 5.0)
    calibration = E.build_count_calibration(pairs(D, M))
    t.close(calibration["agreement"]["rmse"], 2.0, "unhelpful correction: raw RMSE")
    t.close(calibration["loo"][origin]["rmse"], rmse, "unhelpful correction: leave-one-out RMSE")
    t.close(rmse, 2.543545, "unhelpful correction: leave-one-out RMSE as a number", 1e-6)
    t.equal(calibration["recommended"][origin], False, "a correction with a larger leave-one-out error is not recommended")
    t.check(calibration["warnings"][0].startswith("Correction NOT recommended"), "and that is the first warning")
    rows = E.correct_fit_rows(calibration, "run")
    t.equal(rows[0]["Correction_Recommended"], 0, "the fit table says so")
    t.check("NOT recommended" in rows[0]["Recommendation"], "the fit table carries the recommendation text")
    # Raw counts that already agree perfectly: nothing to gain.
    calibration = E.build_count_calibration(pairs([10, 20, 30, 40], [10, 20, 30, 40]))
    t.equal(calibration["recommended"][origin], False, "perfect raw agreement: correction is not recommended")

    # Warnings about the amount and origin of the calibration data.
    def warned(calibration, words):
        return len([text for text in calibration["warnings"] if words in text]) > 0

    three = E.build_count_calibration(pairs([10, 20, 30], [11, 23, 32]))
    t.check(warned(three, "Only 3 validated fields: the calibration is barely determined"), "three fields: strong warning")
    six = E.build_count_calibration(pairs([10, 20, 30, 40, 50, 60], [11, 23, 32, 45, 54, 67]))
    t.check(warned(six, "Only 6 validated fields; 8 or more"), "six fields: mild warning")
    t.check(not warned(six, "barely determined"), "six fields: no strong warning")
    eight = E.build_count_calibration(pairs([10, 20, 30, 40, 50, 60, 70, 80], [11, 23, 32, 45, 54, 67, 76, 89]))
    t.check(not warned(eight, "validated fields"), "eight fields: no warning about their number")
    t.check(warned(eight, "number of images behind the calibration fields is not known"), "plain counts: the number of images is unknown")
    t.check(warned(eight, "settings of the calibration fields are not recorded"), "plain counts: the settings are unknown")
    tagged = pairs([10, 20, 30, 40], [11, 23, 32, 45])
    for index in range(4):
        tagged[index]["source_image"] = "one.tif"
        tagged[index]["parent_roi"] = "field" + str(index // 2)
        tagged[index]["frame"] = "1"
    calibration = E.build_count_calibration(tagged)
    t.check(warned(calibration, "come from 1 image(s)"), "fields of a single image are warned about")
    t.check(warned(calibration, "2 calibration pair(s) repeat a field"), "a field entered twice is warned about")
    t.check(warned(E.build_count_calibration(pairs(CORRECT_TEST_D, CORRECT_TEST_M), line), "cannot be told apart from zero"),
            "line model: an intercept whose interval holds zero is pointed out")
    # Intercept clearly not zero: M = 20 + D with a little scatter.
    shifted = E.build_count_calibration(pairs([10, 20, 30, 40, 50, 60], [30, 41, 50, 59, 70, 81]))
    t.check(warned(shifted, "intercept of the line model differs from zero"), "proportional model: a non-zero intercept is pointed out")

    # The calibration ID names the data and the model, not the order of the pairs.
    first = E.build_count_calibration(pairs([10, 20, 30], [11, 23, 32]))
    second = E.build_count_calibration(pairs([30, 10, 20], [32, 11, 23]))
    third = E.build_count_calibration(pairs([10, 20, 30], [11, 23, 32]), line)
    t.equal(first["id"], second["id"], "the calibration ID does not depend on the order of the fields")
    t.check(first["id"] != third["id"], "another model gives another calibration ID")
    t.check(first["id"].startswith("CAL_") and len(first["id"]) == 16, "the calibration ID is short")


# -----------------------------------------------------------------------------
# 4. Reading tables
# -----------------------------------------------------------------------------

def test_correct_reads_csv(E, t, tmp):
    folder = E.unique_dir(tmp, "correct_csv")
    tricky = [
        "caf" + CORRECT_TEST_EACUTE + " " + MICRO + "m", 'say "hi", twice', "two\nlines", "", "semi;colon", " padded "
    ]
    path = os.path.join(folder, "tricky.csv")
    E.write_csv(path, ["Name", "Value", "Empty"], [{"Name": name, "Value": 1.5} for name in tricky])
    table = E.correct_read_csv(path)
    t.equal(table["headers"], ["Name", "Value", "Empty"], "the header written by write_csv is read back")
    t.equal(correct_test_column(table, "Name"), tricky, "non-ASCII text, quotes, commas and line breaks survive")
    t.equal(correct_test_column(table, "Value"), ["1.5"] * len(tricky), "numbers are read as text")
    t.equal(table["delimiter"], ",", "write_csv files are comma separated")

    # A plain table typed by a user: no quotes, Windows line ends, a
    # byte-order mark, a blank line, a short row.
    path = os.path.join(folder, "plain.csv")
    correct_test_write_bytes(path, b"\xef\xbb\xbfDetected_Count,Manual_Count,Note\r\n10,12,a\r\n\r\n20,25\r\n30,31,c")
    table = E.correct_read_csv(path)
    t.equal(table["headers"], ["Detected_Count", "Manual_Count", "Note"], "a byte-order mark is not part of the first header")
    t.equal(table["rows"], [["10", "12", "a"], ["20", "25", ""], ["30", "31", "c"]], "blank lines are dropped and short rows padded")
    collected = E.collect_calibration_pairs(csv_path=path)
    t.equal([(pair["detected"], pair["manual"]) for pair in collected["pairs"]], [(10.0, 12.0), (20.0, 25.0), (30.0, 31.0)],
            "a simple table gives the calibration pairs")
    t.equal((collected["channel"], collected["channel_name"]), ("", ""), "a table without channel columns names no channel")

    # Semicolons and decimal commas, as spreadsheets write in some regions.
    path = os.path.join(folder, "semicolon.csv")
    correct_test_write_bytes(path, b"Parent_ROI_Name;Detected_Count;Manual_Count\nf1;10,0;12,5\nf2;20;25\nf3;;31\nf4;40;46\n")
    table = E.correct_read_csv(path)
    t.equal(table["delimiter"], ";", "a semicolon-separated table is recognised")
    collected = E.collect_calibration_pairs(csv_path=path)
    t.equal([(pair["detected"], pair["manual"]) for pair in collected["pairs"]], [(10.0, 12.5), (20.0, 25.0), (40.0, 46.0)],
            "decimal commas are read; a row without a count is skipped")
    t.equal(len(collected["skipped"]), 1, "the skipped row is counted")
    t.check(len([note for note in collected["notes"] if "skipped" in note]) == 1, "and reported")

    path = os.path.join(folder, "too_long.csv")
    correct_test_write_bytes(path, b"Detected_Count,Manual_Count\n10,12,99\n")
    t.check(correct_test_raises(lambda: E.correct_read_csv(path)), "a row with more fields than the header is refused")
    path = os.path.join(folder, "no_columns.csv")
    correct_test_write_bytes(path, b"A,B\n1,2\n3,4\n5,6\n")
    t.check(correct_test_raises(lambda: E.collect_calibration_pairs(csv_path=path)), "a table without the two count columns is refused")
    t.check(correct_test_raises(lambda: E.collect_calibration_pairs(csv_path=os.path.join(folder, "missing.csv"))), "a missing table is refused")
    t.check(correct_test_raises(lambda: E.collect_calibration_pairs(folder=E.unique_dir(folder, "empty"))), "a folder without validation runs is refused")

    t.equal([E.correct_channel_number(value) for value in [2, 2.0, "2", "C2", "c2", "C2_DAPI", "DAPI", ""]],
            ["2", "2", "2", "2", "2", "2", "", ""], "channel numbers are understood in every spelling")
    headers = ["ROI_Area", "C1_DAPI_Object_Count", "C3_NeuN_Object_Count", "C1_DAPI_Edge_Touching_Object_Count", "C1_DAPI_Excluded_Edge_Object_Count"]
    t.equal(E.correct_find_count_column(headers, "3"), "C3_NeuN_Object_Count", "the count column is found by channel number")
    t.equal(E.correct_find_count_column(headers, "", None, "DAPI"), "C1_DAPI_Object_Count", "or by channel name")
    t.check(correct_test_raises(lambda: E.correct_find_count_column(headers, "")), "two candidate columns and no channel: refused")
    t.check(correct_test_raises(lambda: E.correct_find_count_column(headers, "2")), "no column of the channel: refused")
    t.equal(E.corrected_column_names("C1_DAPI_Object_Count")["order"], [
        "C1_DAPI_Object_Count_Corrected", "C1_DAPI_Object_Count_Corrected_Lower95", "C1_DAPI_Object_Count_Corrected_Upper95",
        "C1_DAPI_Object_Density_Per_Area_Corrected", "C1_DAPI_Object_Density_Per_Area_Corrected_Lower95",
        "C1_DAPI_Object_Density_Per_Area_Corrected_Upper95",
        "C1_DAPI_Object_Count_Correction_Extrapolated", "C1_DAPI_Object_Count_Correction_Clipped", "Calibration_ID",
        "Correction_Run_ID", "Correction_Application_Version"
    ], "the appended columns and their order")
    # A further count column of the same channel (as another module might
    # add) is not confused with the engine's own column when the channel name
    # is known, and is refused as ambiguous when it is not.
    headers = ["C1_DAPI_Object_Count", "C1_DAPI_Deep_Object_Count", "C2_NeuN_Object_Count"]
    t.equal(E.correct_find_count_column(headers, "1", None, "DAPI"), "C1_DAPI_Object_Count", "number and name give the engine's own count column")
    t.check(correct_test_raises(lambda: E.correct_find_count_column(headers, "1")), "two count columns of one channel and no name: refused")
    t.equal(E.correct_find_count_column(headers, "1", "C1_DAPI_Deep_Object_Count"), "C1_DAPI_Deep_Object_Count", "a column named outright is used")
    t.check(correct_test_raises(lambda: E.correct_find_count_column(headers, "1", "No_Such_Column")), "a named column that does not exist is refused")


# -----------------------------------------------------------------------------
# 5. From real validation runs to a corrected roi_summary.csv
# -----------------------------------------------------------------------------

def correct_test_field(singles, pairs):
    """Disks and manual marks of a synthetic field. A single nucleus is one
    disk with one mark. A merged pair is two disks of radius 6 whose centres
    are 8 pixels apart: they overlap, so the detection finds ONE object, but
    the annotator marks both. A field with s singles and p pairs therefore
    has s + p detected objects and s + 2p manual marks."""
    disks = []
    for position in singles:
        disks.append(position)
    for position in pairs:
        disks.append(position)
        disks.append((position[0] + 8, position[1]))
    return disks


def correct_test_l_roi(name):
    """An L-shaped polygon: the top band (y < 95) of the 240 x 200 image plus
    the left half below it. Not a rectangle, and it touches three borders."""
    polygon = Polygon()
    for point in [(0, 0), (240, 0), (240, 95), (120, 95), (120, 200), (0, 200)]:
        polygon.addPoint(point[0], point[1])
    roi = PolygonRoi(polygon, Roi.POLYGON)
    roi.setName(name)
    return roi


def correct_test_named(roi, name):
    roi.setName(name)
    return roi


def correct_test_validation(E, root, title, disks, rois, config, with_settings=True):
    """One validation run, written the way the window writes it."""
    imp = make_image(title, slices=1, disks_per_frame=[disks])
    engine = E.AnalysisEngine()
    prepared = engine.prepare_single_channel(imp, config, 1)
    run_dir = E.unique_dir(root, E.sanitize_name(E.strip_extension(title)) + "_CellQuant_Validation_test")
    try:
        if with_settings:
            E.write_settings_properties(os.path.join(run_dir, "settings.properties"), E.default_global_options(), [config])
        # one mark on the centre of every disk
        return E.validate_objects_against_points(engine, prepared, config, rois, list(disks), run_dir, title, 1)
    finally:
        prepared.close()
        imp.close()


def correct_test_calibration_runs(E, root, with_settings=True):
    """Three validation runs with two parent ROIs each. By construction
    (singles s, merged pairs p -> detected s + p, manual s + 2p):
        run a: left  s2 p2 -> D 4, M 6      right s4 p1 -> D 5, M 6
        run b: left  s1 p3 -> D 4, M 7      right s3 p3 -> D 6, M 9
        run c: L     s6 p1 -> D 7, M 8      box   s2 p4 -> D 6, M 10
    Run c uses the L-shaped polygon, and one of its singles is a disk cut by
    the right image border (centre x = 237 in a 240 pixel wide image)."""
    config = nuclei_config(E)
    left = lambda: correct_test_named(Roi(0, 0, 120, HEIGHT), "left")
    right = lambda: correct_test_named(Roi(120, 0, 120, HEIGHT), "right")
    correct_test_validation(E, root, "cal_a.tif", correct_test_field(
        [(20, 20), (55, 20), (140, 20), (175, 20), (210, 20), (140, 50)],
        [(20, 50), (55, 50), (175, 50)]), [left(), right()], config, with_settings)
    correct_test_validation(E, root, "cal_b.tif", correct_test_field(
        [(20, 20), (140, 20), (175, 20), (210, 20)],
        [(55, 20), (20, 50), (55, 50), (140, 50), (175, 50), (140, 80)]), [left(), right()], config, with_settings)
    return correct_test_validation(E, root, "cal_c " + MICRO + ".tif", correct_test_field(
        [(20, 20), (55, 20), (90, 20), (140, 20), (20, 110), (237, 50), (140, 110), (175, 110)],
        [(55, 110), (210, 110), (140, 140), (175, 140), (210, 140)]),
        [correct_test_l_roi("L " + CORRECT_TEST_EACUTE), correct_test_named(Roi(120, 95, 120, 105), "box")], config, with_settings)


def correct_test_table_run(E, out_root, overrides=None):
    """The image whose counts are corrected. Detected objects by construction:
    whole image 7 (one of them cut by the right image border), L-shaped
    polygon 6, bottom-right box 1 (the border object), tiny empty ROI 0."""
    disks = correct_test_field([(20, 20), (55, 50), (140, 20), (20, 140), (237, 150)], [(90, 80), (55, 170)])
    imp = make_image("section " + MICRO + ".tif", slices=1, disks_per_frame=[disks])
    rois = [
        correct_test_named(Roi(0, 0, WIDTH, HEIGHT), "whole"), correct_test_l_roi("L " + CORRECT_TEST_EACUTE),
        correct_test_named(Roi(120, 95, 120, 105), "box"), correct_test_named(Roi(2, 2, 3, 3), "empty")
    ]
    try:
        return run_image(E, imp, [nuclei_config(E, overrides), marker_config(E)], rois, out_root)
    finally:
        imp.close()


def test_correct_from_validation_runs(E, t, tmp):
    root = E.unique_dir(tmp, "correct_runs")
    cal_c1 = os.path.join(root, "calibration", "dapi")
    cal_c2 = os.path.join(root, "calibration", "map2")
    os.makedirs(cal_c1)
    os.makedirs(cal_c2)
    last = correct_test_calibration_runs(E, cal_c1)
    t.equal([(row["Detected_Count"], row["Manual_Count"]) for row in last["rows"]], [(7, 8), (6, 10), (13, 18)],
            "the validation of run c gives the counts built into the image, and a pooled row")

    # --- collecting the pairs ---
    collected = E.collect_calibration_pairs(cal_c1)
    pairs = collected["pairs"]
    t.equal(sorted([(pair["detected"], pair["manual"]) for pair in pairs]),
            [(4, 6), (4, 7), (5, 6), (6, 9), (6, 10), (7, 8)], "one pair per validated parent ROI; pooled rows are skipped")
    t.equal(len([pair for pair in pairs if pair["parent_roi"] == "ALL_ROIS_POOLED"]), 0, "no pooled row among the pairs")
    t.equal((collected["channel"], collected["channel_name"], collected["files"]), ("1", "DAPI", 3), "channel and file count are reported")
    t.check("L " + CORRECT_TEST_EACUTE in [pair["parent_roi"] for pair in pairs], "a non-ASCII ROI name is read")
    t.check("cal_c " + MICRO + ".tif" in [pair["source_image"] for pair in pairs], "a non-ASCII image name is read")
    t.equal(len(set([pair["fp_full"] for pair in pairs])), 1, "runs with the same settings share one settings fingerprint")
    t.check(len(pairs[0]["fp_full"]) == 12 and len(pairs[0]["fp_partial"]) == 12, "both fingerprints are recorded")
    t.check(correct_test_raises(lambda: E.collect_calibration_pairs(os.path.join(cal_c1, os.listdir(cal_c1)[0]) + "_missing")),
            "a calibration folder that does not exist is refused")

    # A second channel in the same parent folder: C2 as objects (the marker
    # rectangle is one object, marked three times -> D 1, M 3).
    marker = marker_config(E, {"analysis_mode": "Objects only"})
    imp = make_image("cal_marker.tif", slices=1)
    engine = E.AnalysisEngine()
    prepared = engine.prepare_single_channel(imp, marker, 1)
    E.validate_objects_against_points(engine, prepared, marker, [E.whole_image_roi(imp)], [(50, 60), (90, 90), (120, 120)],
                                      E.unique_dir(cal_c2, "marker_CellQuant_Validation_test"), "cal_marker.tif", 1)
    prepared.close()
    imp.close()
    both = os.path.join(root, "calibration")
    t.check(correct_test_raises(lambda: E.collect_calibration_pairs(both)), "two channels and no choice: refused")
    for wanted in ["1", 1, "C1", "DAPI", "C1_DAPI", "c1_dapi"]:
        t.equal(len(E.collect_calibration_pairs(both, None, wanted)["pairs"]), 6, "channel chosen as " + repr(wanted))
    t.equal(len(E.collect_calibration_pairs(both, None, "MAP2")["pairs"]), 1, "the other channel can be chosen too")
    t.check(correct_test_raises(lambda: E.run_count_correction(calibration_folder=cal_c2, output_root=os.path.join(root, "refused"))),
            "a single validated field is refused")
    t.check(not os.path.exists(os.path.join(root, "refused")), "a refused calibration writes nothing")

    # --- the table to correct ---
    result = correct_test_table_run(E, os.path.join(root, "runs"))
    t.equal([row["C1_DAPI_Object_Count"] for row in result["summary_rows"]], [7, 6, 1, 0],
            "the image to correct holds the object counts built into it")
    t.equal(result["summary_rows"][0]["C1_DAPI_Edge_Touching_Object_Count"], 1, "one object touches the image border")
    t.close(result["summary_rows"][0]["ROI_Area"], WIDTH * HEIGHT * PIXEL * PIXEL, "whole-image area")
    summary_path = os.path.join(result["run_dir"], "roi_summary.csv")
    before_bytes = correct_test_read_bytes(summary_path)
    before_files = sorted(os.listdir(result["run_dir"]))

    done = E.run_count_correction(calibration_folder=cal_c1, table_path=summary_path, output_root=os.path.join(root, "corrected"))
    t.equal(correct_test_read_bytes(summary_path), before_bytes, "the original roi_summary.csv is untouched")
    t.equal(sorted(os.listdir(result["run_dir"])), before_files, "nothing is added to the run folder")
    out_dir = done["out_dir"]
    t.equal(sorted(os.listdir(out_dir)), ["count_calibration.csv", "count_calibration_fit.csv", "roi_summary_corrected.csv"],
            "three new files are written")
    t.check("_CellQuant_" in os.path.basename(out_dir), "the output folder is recognised as an application folder by batches")

    # Proportional model through the six pairs, by hand:
    #   sum D D = 16 + 25 + 16 + 36 + 49 + 36 = 178
    #   sum D M = 24 + 30 + 28 + 54 + 56 + 60 = 252
    #   sum M M = 36 + 36 + 49 + 81 + 64 + 100 = 366
    #   slope = 252 / 178 = 1.415730
    #   SSE = 366 - 252^2 / 178 = 9.235955 on 5 degrees of freedom, SD = 1.359114
    slope = 252.0 / 178.0
    sd = math.sqrt((366.0 - 252.0 * 252.0 / 178.0) / 5.0)
    fit = done["calibration"]["fit"]
    t.close(fit["slope"], slope, "slope fitted from the validation files")
    t.close(fit["residual_sd"], sd, "residual SD fitted from the validation files")
    t.close(fit["residual_sd"], 1.359114, "residual SD as a number", 1e-6)
    t.equal((fit["n"], fit["df"]), (6, 5), "six fields, five degrees of freedom")
    t.close(fit["t"], CORRECT_TEST_T5, "t quantile on 5 degrees of freedom")
    # raw RMSE = sqrt((4 + 1 + 9 + 9 + 1 + 16) / 6) = 2.582. Leave-one-out:
    # the slope without pair (D, M) is (252 - D M) / (178 - D D), which gives
    # the errors -0.370 1.255 -1.469 -0.634 2.636 -1.887 and an RMSE of 1.570.
    # That is smaller than 2.582, so the correction is recommended.
    built = [(4, 6), (5, 6), (4, 7), (6, 9), (7, 8), (6, 10)]
    errors = [(252.0 - d * m) / (178.0 - d * d) * d - m for d, m in built]
    loo_rmse = math.sqrt(sum([error * error for error in errors]) / 6.0)
    t.close(loo_rmse, 1.570, "leave-one-out RMSE worked by hand", 1e-3)
    t.close(done["calibration"]["agreement"]["rmse"], math.sqrt(40.0 / 6.0), "raw RMSE of the validation fields")
    t.close(done["calibration"]["loo"][E.CORRECT_MODEL_ORIGIN]["rmse"], loo_rmse, "leave-one-out RMSE of the corrected counts")
    t.equal(done["calibration"]["recommended"][E.CORRECT_MODEL_ORIGIN], True, "the correction is recommended")
    t.equal(done["calibration"]["n_images"], 3, "the fields come from three images")

    # --- the corrected table ---
    original_lines = read_utf8(summary_path).split("\n")
    corrected_lines = read_utf8(done["corrected_path"]).split("\n")
    names = E.corrected_column_names("C1_DAPI_Object_Count")
    t.equal(len(corrected_lines), len(original_lines), "the corrected table has the rows of the original")
    t.equal(corrected_lines[0], original_lines[0] + "," + ",".join(['"' + name + '"' for name in names["order"]]),
            "the header is the original header with the new columns appended")
    for index in range(1, 5):
        t.check(corrected_lines[index].startswith(original_lines[index] + ","),
                "row " + str(index) + ": every original column is unchanged and in place")
    table = E.correct_read_csv(done["corrected_path"])
    t.equal(correct_test_column(table, "Parent_ROI_Name"), ["whole", "L " + CORRECT_TEST_EACUTE, "box", "empty"], "ROI names survive, non-ASCII included")
    t.equal(correct_test_column(table, "Source_Image"), ["section " + MICRO + ".tif"] * 4, "the non-ASCII image name survives")
    t.equal(set(correct_test_column(table, "Calibration_ID")), set([done["calibration"]["id"]]), "every row names the calibration")
    t.equal(set(correct_test_column(table, "Correction_Run_ID")), set([os.path.basename(out_dir)]), "every row names the correction folder")
    t.equal(set(correct_test_column(table, "Correction_Application_Version")), set([E.APP_VERSION]), "and the version that corrected it")
    t.equal(set(correct_test_column(table, "Run_ID")), set([os.path.basename(result["run_dir"])]), "the Run_ID of the analysis is kept")
    count = [float(value) for value in correct_test_column(table, names["count"])]
    lower = [float(value) for value in correct_test_column(table, names["count_low"])]
    upper = [float(value) for value in correct_test_column(table, names["count_high"])]
    # corrected = slope x detected: 7 -> 9.910112, 6 -> 8.494382, 1 -> 1.415730, 0 -> 0
    for index in range(4):
        t.close(count[index], slope * [7, 6, 1, 0][index], "corrected count of row " + str(index + 1))
    t.close(count[1], 8.494382, "corrected count of the polygon ROI as a number", 1e-6)
    # prediction interval: +/- t x SD x sqrt(1 + D^2 / 178)
    #   t x SD = 2.570582 x 1.359114 = 3.493714
    #   whole (D 7): 3.493714 x sqrt(1 + 49/178) = 3.493714 x 1.129284 = 3.945393
    #   L     (D 6): 3.493714 x sqrt(1 + 36/178) = 3.493714 x 1.096470 = 3.830754
    half = [CORRECT_TEST_T5 * sd * math.sqrt(1.0 + d * d / 178.0) for d in [7.0, 6.0, 1.0, 0.0]]
    t.close(half[1], 3.830754, "half width for the polygon ROI as a number", 1e-6)
    for index in range(2):
        t.close(lower[index], count[index] - half[index], "lower prediction limit of row " + str(index + 1))
        t.close(upper[index], count[index] + half[index], "upper prediction limit of row " + str(index + 1))
    # box (D 1): 1.4157 - 3.5035 is below zero -> lower limit 0; empty (D 0): 0 +/- 3.4937
    t.equal((lower[2], lower[3]), (0.0, 0.0), "lower limits are not negative")
    t.close(upper[2], slope + half[2], "upper limit of the one-object ROI")
    t.close(upper[3], CORRECT_TEST_T5 * sd, "an empty ROI still has an upper limit: t x SD")
    t.equal(correct_test_column(table, names["extrapolated"]), ["0", "0", "1", "1"],
            "counts of 1 and 0 are below the calibrated range 4 to 7 and flagged; 7 and 6 are not")
    t.equal(correct_test_column(table, names["clipped"]), ["0", "0", "0", "0"], "nothing was below zero")
    # density = corrected count / ROI_Area; the whole image is 240 x 200 x 0.25 = 12000 square micrometers
    area = [float(value) for value in correct_test_column(table, "ROI_Area")]
    t.close(float(correct_test_column(table, names["density"])[0]), 7.0 * slope / 12000.0, "corrected density of the whole image")
    for index in range(4):
        t.close(float(correct_test_column(table, names["density"])[index]), count[index] / area[index], "density of row " + str(index + 1))
        t.close(float(correct_test_column(table, names["density_low"])[index]), lower[index] / area[index], "density lower limit of row " + str(index + 1))
        t.close(float(correct_test_column(table, names["density_high"])[index]), upper[index] / area[index], "density upper limit of row " + str(index + 1))
    t.equal((done["corrected"]["corrected"], done["corrected"]["extrapolated"], done["corrected"]["clipped"]), (4, 2, 0), "flag totals")
    t.equal(done["corrected"]["settings_match"], "SAME", "the table was made with the settings of the validation runs")
    t.check(len([text for text in done["warnings"] if "OUTSIDE the calibrated range" in text]) == 1, "extrapolation is warned about")
    t.check(len([text for text in done["warnings"] if "DIFFERENT detection settings" in text]) == 0, "no settings warning")

    # --- the two calibration tables ---
    pairs_table = E.correct_read_csv(done["pairs_path"])
    t.equal(pairs_table["headers"], E.COUNT_CALIBRATION_HEADERS, "count_calibration.csv has its fixed header")
    t.equal(len(pairs_table["rows"]), 6, "one row per calibration pair")
    t.equal(set(correct_test_column(pairs_table, "Run_ID")), set([os.path.basename(out_dir)]), "Run_ID is the output folder")
    t.equal(set(correct_test_column(pairs_table, "Application_Version")), set([E.APP_VERSION]), "the application version is recorded")
    t.equal(set(correct_test_column(pairs_table, "Validation_Application_Version")), set([E.APP_VERSION]), "and that of the validation runs")
    t.check("L " + CORRECT_TEST_EACUTE in correct_test_column(pairs_table, "Parent_ROI_Name"), "non-ASCII ROI names reach count_calibration.csv")
    t.equal(sorted(zip(correct_test_column(pairs_table, "Detected_Count"), correct_test_column(pairs_table, "Manual_Count"))),
            [("4", "6"), ("4", "7"), ("5", "6"), ("6", "10"), ("6", "9"), ("7", "8")], "whole counts are written as whole numbers")
    for row_index in range(6):
        d = float(correct_test_column(pairs_table, "Detected_Count")[row_index])
        m = float(correct_test_column(pairs_table, "Manual_Count")[row_index])
        t.close(float(correct_test_column(pairs_table, "Fitted_Manual_Count")[row_index]), slope * d, "fitted value of pair " + str(row_index + 1))
        t.close(float(correct_test_column(pairs_table, "Residual_Manual_Minus_Fitted")[row_index]), m - slope * d, "residual of pair " + str(row_index + 1))
        # leave-one-out: slope without this pair is (252 - D M) / (178 - D D)
        t.close(float(correct_test_column(pairs_table, "LeaveOneOut_Corrected_Count")[row_index]), (252.0 - d * m) / (178.0 - d * d) * d,
                "leave-one-out corrected count of pair " + str(row_index + 1))
    fit_table = E.correct_read_csv(done["fit_path"])
    t.equal(fit_table["headers"], E.COUNT_CALIBRATION_FIT_HEADERS, "count_calibration_fit.csv has its fixed header")
    t.equal(correct_test_column(fit_table, "Model"), E.CORRECT_MODELS, "one row per model")
    t.equal(correct_test_column(fit_table, "Selected"), ["1", "0"], "the applied model is marked")
    t.equal(correct_test_column(fit_table, "N_Fields"), ["6", "6"], "number of fields")
    t.equal(correct_test_column(fit_table, "N_Images"), ["3", "3"], "number of images")
    t.close(float(correct_test_column(fit_table, "Slope")[0]), slope, "slope in the fit table")
    t.equal(correct_test_column(fit_table, "Settings_Match"), ["SAME", "SAME"], "settings verdict in the fit table")
    t.equal(correct_test_column(fit_table, "Rows_Extrapolated")[0], "2", "extrapolated rows in the fit table")
    t.equal(correct_test_column(fit_table, "Correction_Recommended")[0], "1", "recommendation in the fit table")
    t.check("counting error only" in correct_test_column(fit_table, "Interval_Note")[0], "the scope of the interval is stated in the table")
    t.check(correct_test_column(fit_table, "Intercept_Lower95")[1] != "" and correct_test_column(fit_table, "Intercept_Lower95")[0] == "",
            "the intercept interval is reported for the line model only")

    # --- the same table made with other settings is recognised ---
    other = correct_test_table_run(E, os.path.join(root, "runs_other"), {"min_size": 5})
    t.equal([row["C1_DAPI_Object_Count"] for row in other["summary_rows"]], [7, 6, 1, 0], "the counts happen to be the same")
    changed = E.correct_count_table(os.path.join(other["run_dir"], "roi_summary.csv"), done["calibration"], os.path.join(root, "corrected_other"))
    t.equal(changed["settings_match"], "DIFFERENT", "a table made with another minimum size is recognised")
    t.check(len([text for text in changed["warnings"] if "DIFFERENT detection settings" in text]) == 1, "and warned about")
    t.check(correct_test_raises(lambda: E.correct_count_table(os.path.join(other["run_dir"], "roi_summary.csv"), done["calibration"], os.path.join(root, "corrected_other"))),
            "an existing corrected table is not overwritten")
    t.check(correct_test_raises(lambda: E.correct_count_table(changed["path"], done["calibration"], os.path.join(root, "corrected_twice"))),
            "a corrected table is not corrected a second time")

    # --- validation runs without settings.properties: thresholds.csv is used ---
    bare = os.path.join(root, "calibration_bare")
    os.makedirs(bare)
    correct_test_calibration_runs(E, bare, False)
    partial = E.build_count_calibration(E.collect_calibration_pairs(bare)["pairs"], channel="1", channel_name="DAPI")
    t.check("thresholds.csv" in partial["fingerprint_source"], "without settings.properties the fingerprint comes from thresholds.csv")
    same = E.correct_count_table(summary_path, partial, os.path.join(root, "corrected_partial"))
    t.equal(same["settings_match"], "SAME_SEGMENTATION_SETTINGS_ONLY", "segmentation settings of thresholds.csv agree with the run's settings.properties")
    differs = E.correct_count_table(os.path.join(other["run_dir"], "roi_summary.csv"), partial, os.path.join(root, "corrected_partial_other"))
    t.equal(differs["settings_match"], "DIFFERENT", "and a different minimum size is still recognised")
    t.equal(partial["id"], done["calibration"]["id"], "the same pairs, model and channel give the same calibration ID")


def test_correct_changes_nothing_else(E, t, tmp):
    """The correction is a separate step on tables. Loading the module must
    leave every analysis table as it was, and profiles of older versions
    must keep loading."""
    t.equal(E.ENGINE_MODULE_VERSIONS.get("cellquantcorrect.py"), E.APP_VERSION, "the module is loaded and carries the application version")
    options = E.default_global_options()
    t.equal(options["correct_model"], E.CORRECT_MODEL_ORIGIN, "default model: proportional")
    t.equal(options["correct_weighting"], E.CORRECT_WEIGHT_NONE, "default weighting: ordinary least squares")
    configs = [nuclei_config(E), marker_config(E)]
    engine = E.AnalysisEngine()
    t.equal(engine.build_summary_headers(configs, options), engine.build_summary_headers(configs), "default settings add no summary column")
    t.equal(engine.build_object_headers(configs, options), engine.build_object_headers(configs), "default settings add no object column")
    for model in E.CORRECT_MODELS:
        changed = dict(options)
        changed["correct_model"] = model
        changed["correct_weighting"] = E.CORRECT_WEIGHT_COUNT
        t.equal(engine.build_summary_headers(configs, changed), engine.build_summary_headers(configs),
                "the correction settings never change the analysis tables (" + model + ")")
    # A profile saved before this module existed has no correct.* keys.
    path = os.path.join(tmp, "correct_profile.properties")
    options["correct_model"] = E.CORRECT_MODEL_LINE
    E.write_settings_properties(path, options, configs)
    props = E.read_properties_file(path)
    t.equal(E.global_options_from_properties(props)["correct_model"], E.CORRECT_MODEL_LINE, "the model is saved in a profile")
    props.remove("correct.model")
    props.remove("correct.weighting")
    loaded = E.global_options_from_properties(props)
    t.equal((loaded["correct_model"], loaded["correct_weighting"]), (E.CORRECT_MODEL_ORIGIN, E.CORRECT_WEIGHT_NONE),
            "a profile without the new keys loads with the defaults")
    # run_count_correction takes its model from the options when none is given.
    folder = E.unique_dir(tmp, "correct_options")
    table = os.path.join(folder, "pairs.csv")
    correct_test_write_bytes(table, b"Detected_Count,Manual_Count\n10,12\n20,25\n30,31\n40,46\n50,56\n")
    done = E.run_count_correction(calibration_csv=table, global_options={"correct_model": E.CORRECT_MODEL_LINE}, output_root=folder)
    t.equal(done["calibration"]["model"], E.CORRECT_MODEL_LINE, "the model of the profile is used")
    t.close(done["calibration"]["fit"]["slope"], 1.09, "and fitted (worked example, line model)")
    t.equal(done["corrected"], None, "without a table only the calibration is written")
    t.check(os.path.isfile(done["fit_path"]) and os.path.isfile(done["pairs_path"]), "both calibration tables are written")
    # A plain table of counts without ROI_Area: corrected, density left empty.
    counts = os.path.join(folder, "counts.csv")
    correct_test_write_bytes(counts, b"Field,C1_DAPI_Object_Count\na,30\nb,\nc,500\n")
    result = E.correct_count_table(counts, E.build_count_calibration(E.correct_pairs_from_counts(CORRECT_TEST_D, CORRECT_TEST_M)), E.unique_dir(folder, "out"))
    rows = result["rows"]
    t.close(rows[0]["C1_DAPI_Object_Count_Corrected"], 30.0 * 6190.0 / 5500.0, "a plain table is corrected")
    t.check("C1_DAPI_Object_Density_Per_Area_Corrected" not in rows[0], "without ROI_Area there is no density")
    t.check("C1_DAPI_Object_Count_Corrected" not in rows[1], "a row without a count stays uncorrected")
    t.equal((result["corrected"], result["blank"], result["extrapolated"]), (2, 1, 1), "row totals")
    t.equal(rows[2]["C1_DAPI_Object_Count_Correction_Extrapolated"], 1, "500 detected is far above the calibrated range 10 to 50")
    t.check(len([text for text in result["warnings"] if "more objects than any calibration field" in text]) == 1,
            "ROIs larger than the calibration fields get their own explanation")


def correct_test_summary_file(E, folder, counts, config=None):
    """A validation folder written by hand: validation_summary.csv with one
    row per (detected, manual) pair, each from its own image, and optionally
    the settings.properties of a channel configuration."""
    rows = []
    for index in range(len(counts)):
        rows.append({
            "Run_ID": os.path.basename(folder), "Application_Version": E.APP_VERSION,
            "Source_Image": "image" + str(index) + ".tif", "Frame": 1, "Channel": 1, "Channel_Name": "DAPI",
            "Parent_ROI_Index": 1, "Parent_ROI_Name": "field", "Detected_Count": counts[index][0], "Manual_Count": counts[index][1]
        })
    E.write_csv(os.path.join(folder, "validation_summary.csv"), E.VALIDATION_SUMMARY_HEADERS, rows)
    if config is not None:
        E.write_settings_properties(os.path.join(folder, "settings.properties"), E.default_global_options(), [config])


def test_correct_settings_checks(E, t, tmp):
    """What the settings comparison can and cannot vouch for."""
    root = E.unique_dir(tmp, "correct_settings")
    counts = [(10, 12), (20, 25), (30, 31), (40, 46), (50, 56)]

    def warned(texts, words):
        return len([text for text in texts if words in text]) > 0

    def table_folder(name, config):
        folder = E.unique_dir(root, name)
        correct_test_write_bytes(os.path.join(folder, "roi_summary.csv"), b"ROI_Area,C1_DAPI_Object_Count\n100,30\n")
        if config is not None:
            E.write_settings_properties(os.path.join(folder, "settings.properties"), E.default_global_options(), [config])
        return os.path.join(folder, "roi_summary.csv")

    # A threshold set by eye ('Display minimum') is not part of the profile:
    # the profile can be identical while the thresholds differ.
    by_eye = nuclei_config(E, {"threshold_mode": "Display minimum"})
    cal_eye = E.unique_dir(root, "cal_eye")
    correct_test_summary_file(E, cal_eye, counts, by_eye)
    collected = E.collect_calibration_pairs(cal_eye)
    t.equal([pair["by_eye"] for pair in collected["pairs"]], [True] * 5, "fields thresholded by eye are recognised")
    calibration = E.build_count_calibration(collected["pairs"], channel="1", channel_name="DAPI")
    t.check(warned(calibration["warnings"], "set by eye"), "a calibration from by-eye thresholds is warned about")
    done = E.correct_count_table(table_folder("table_eye", by_eye), calibration, E.unique_dir(root, "out_eye"))
    t.equal(done["settings_match"], "THRESHOLD_SET_BY_EYE", "the same profile with a by-eye threshold is not called the same")
    t.check(warned(done["warnings"], "set by eye"), "and the table gets its own warning")
    # worked example of the file header: slope 6190 / 5500, so 30 detected -> 33.763636
    t.close(done["rows"][0]["C1_DAPI_Object_Count_Corrected"], 30.0 * 6190.0 / 5500.0, "the count is corrected all the same")

    # Fixed thresholds: identical settings are SAME, no by-eye warning.
    fixed = nuclei_config(E)
    cal_fixed = E.unique_dir(root, "cal_fixed")
    correct_test_summary_file(E, cal_fixed, counts, fixed)
    collected = E.collect_calibration_pairs(cal_fixed)
    t.equal([pair["by_eye"] for pair in collected["pairs"]], [False] * 5, "a fixed threshold is not flagged")
    calibration = E.build_count_calibration(collected["pairs"], channel="1", channel_name="DAPI")
    t.check(not warned(calibration["warnings"], "set by eye"), "no by-eye warning for a fixed threshold")
    t.check(not warned(calibration["warnings"], "come from"), "five fields from five images: no warning about the number of images")
    same = E.correct_count_table(table_folder("table_same", fixed), calibration, E.unique_dir(root, "out_same"))
    t.equal(same["settings_match"], "SAME", "identical fixed settings are recognised as the same")
    # Display settings do not change detection and must not change the verdict;
    # the manual threshold does.
    look = E.correct_count_table(table_folder("table_lut", nuclei_config(E, {"lut": "Cyan"})), calibration, E.unique_dir(root, "out_lut"))
    t.equal(look["settings_match"], "SAME", "another display colour is still the same detection")
    moved = E.correct_count_table(table_folder("table_moved", nuclei_config(E, {"manual_low": 1200})), calibration, E.unique_dir(root, "out_moved"))
    t.equal(moved["settings_match"], "DIFFERENT", "another manual threshold is a different detection")
    eye_table = E.correct_count_table(table_folder("table_eye2", by_eye), calibration, E.unique_dir(root, "out_eye2"))
    t.equal(eye_table["settings_match"], "DIFFERENT", "fixed against by-eye threshold is different")
    none = E.correct_count_table(table_folder("table_none", None), calibration, E.unique_dir(root, "out_none"))
    t.equal(none["settings_match"], "NOT_CHECKED", "a table without recorded settings cannot be checked")
    t.check(warned(none["warnings"], "could not be checked"), "and that is said")

    # Validation runs of two different profiles in one folder are mixed.
    mixed_root = E.unique_dir(root, "cal_mixed")
    correct_test_summary_file(E, E.unique_dir(mixed_root, "a"), counts[:3], fixed)
    correct_test_summary_file(E, E.unique_dir(mixed_root, "b"), counts[3:], nuclei_config(E, {"manual_low": 1200}))
    mixed = E.build_count_calibration(E.collect_calibration_pairs(mixed_root)["pairs"], channel="1", channel_name="DAPI")
    t.equal(mixed["fingerprint"], "MIXED", "validation runs of two profiles are recognised")
    t.check(warned(mixed["warnings"], "DIFFERENT detection settings"), "and warned about")

    # A table that cannot be corrected leaves nothing behind.
    out_root = os.path.join(root, "out_refused")
    bad = os.path.join(root, "no_count_column.csv")
    correct_test_write_bytes(bad, b"ROI_Area,Something\n100,30\n")
    t.check(correct_test_raises(lambda: E.run_count_correction(calibration_folder=cal_fixed, table_path=bad, output_root=out_root)),
            "a table without an object-count column is refused")
    t.equal(os.listdir(out_root) if os.path.isdir(out_root) else [], [], "and no empty output folder is left behind")


# -----------------------------------------------------------------------------
# 6. Further known answers added by the independent check
# -----------------------------------------------------------------------------

def correct_test_corner():
    """The triangle (0,0) (140,0) (0,140): non-rectangular, on the top and
    left image borders."""
    polygon = Polygon()
    for point in [(0, 0), (140, 0), (0, 140)]:
        polygon.addPoint(point[0], point[1])
    return correct_test_named(PolygonRoi(polygon, Roi.POLYGON), "corner")


def test_correct_weighted_line_and_line_leave_one_out(E, t, tmp):
    """Two cases the tests above do not reach: the LINE model with weights
    1 / detected count, and leave-one-out of the LINE model, worked with
    exact fractions."""
    line = E.CORRECT_MODEL_LINE
    origin = E.CORRECT_MODEL_ORIGIN
    t2 = 4.302652730    # t(0.975) on 2 degrees of freedom, printed tables

    # --- weighted line: D = 5 10 20 40, M = 8 12 25 44, weights w = 1 / D ---
    # sum w = 1/5 + 1/10 + 1/20 + 1/40 = 3/8
    # sum w D = 4 (each term is 1)       -> weighted mean D = 4 / (3/8) = 32/3
    # sum w M = 1.6 + 1.2 + 1.25 + 1.1 = 5.15 -> weighted mean M = 206/15
    # Sxx_w = sum w D^2 - (sum w D)^2 / sum w = 75 - 16 x 8/3 = 97/3
    # Sxy_w = sum w D M - (sum w D)(sum w M) / sum w = 89 - 4 x 5.15 x 8/3 = 511/15
    # slope = (511/15) / (97/3) = 511/485 = 1.053608
    # intercept = 206/15 - (511/485)(32/3) = 242/97 = 2.494845
    # residuals 0.237113 -1.030928 1.432990 -0.639175;
    # weighted SSE = sum r^2 / D = 447/1940 on 2 degrees of freedom
    # residual SD = sqrt(447/3880) = 0.339420
    # SE(slope) = SD / sqrt(Sxx_w) = 0.059692
    # SE(intercept) = SD x sqrt(1/sum w + mean^2 / Sxx_w) = 0.844166
    # A new field with 30 detected: estimate 242/97 + 30 x 511/485 = 34.103093;
    # its variance is SD^2 x (30 + 1/sum w + (30 - 32/3)^2 / Sxx_w)
    #   = SD^2 x (30 + 8/3 + 3364/9 x 3/97) = SD^2 x 4290/97 = SD^2 x 44.226804
    # half width = t2 x SD x sqrt(44.226804) = 9.712186
    fit = E.fit_count_calibration([5, 10, 20, 40], [8, 12, 25, 44], line, E.CORRECT_WEIGHT_COUNT)
    t.close(fit["slope"], 511.0 / 485.0, "weighted line: slope")
    t.close(fit["intercept"], 242.0 / 97.0, "weighted line: intercept")
    t.close(fit["residual_sd"], 0.339420, "weighted line: residual SD", 1e-5)
    t.close(fit["residual_sd"], math.sqrt(447.0 / 3880.0), "weighted line: residual SD exactly")
    t.close(fit["slope_se"], 0.059692, "weighted line: SE of the slope", 1e-4)
    t.close(fit["intercept_se"], 0.844166, "weighted line: SE of the intercept", 1e-5)
    t.equal(fit["df"], 2, "weighted line: n - 2 degrees of freedom")
    t.close(fit["intercept_low"], 242.0 / 97.0 - t2 * 0.844166, "weighted line: intercept lower limit", 1e-5)
    result = E.predict_corrected_count(fit, 30)
    t.close(result["estimate"], 242.0 / 97.0 + 30.0 * 511.0 / 485.0, "weighted line: corrected count of 30 detected")
    t.close(result["upper"] - result["estimate"], 9.712186, "weighted line: half width of the prediction interval", 1e-6)
    t.close(result["estimate"] - result["lower"], t2 * math.sqrt(447.0 / 3880.0) * math.sqrt(4290.0 / 97.0),
            "weighted line: the interval is symmetric")

    # --- line, ordinary least squares: D = 10 20 30 40, M = 12 19 33 41 ---
    # mean D 25, mean M 26.25, Sxx 500, Sxy 505 -> slope 1.01, intercept 1
    # fitted 11.1 21.2 31.3 41.4, residuals 0.9 -2.2 1.7 -0.4, SSE 8.7 on 2 df
    # SD = sqrt(4.35) = 2.085665; SE(intercept) = SD x sqrt(1/4 + 625/500) = 2.554408
    # intercept CI = 1 +/- 4.302653 x 2.554408 = -9.990730 to 11.990730
    D = [10, 20, 30, 40]
    M = [12, 19, 33, 41]
    fit = E.fit_count_calibration(D, M, line)
    t.close(fit["slope"], 1.01, "line: slope")
    t.close(fit["intercept"], 1.0, "line: intercept")
    t.close(fit["residual_sd"], math.sqrt(4.35), "line: residual SD")
    t.close(fit["intercept_low"], -9.990730, "line: intercept lower limit", 1e-6)
    t.close(fit["intercept_high"], 11.990730, "line: intercept upper limit", 1e-6)
    t.check(fit["aicc"] != fit["aicc"], "line with 4 fields: AICc is not defined (n - k - 1 = 0) and is left empty")
    # Leave-one-out of the line, each fold a line through three points:
    #   without (10,12): (20,19) (30,33) (40,41): slope 220/200 = 1.1, intercept 31 - 33 = -2 -> 9
    #   without (20,19): slope 137/140, intercept 18/7 -> 18/7 + 20 x 137/140 = 155/7 = 22.142857
    #   without (30,33): slope 69/70, intercept 1 -> 1 + 30 x 69/70 = 30.571429
    #   without (40,41): (10,12) (20,19) (30,33): slope 210/200 = 1.05, intercept 64/3 - 21 = 1/3 -> 127/3
    # errors -3, 22/7, -17/7, 4/3: RMSE = sqrt((9 + 484/49 + 289/49 + 16/9) / 4) = 2.576494
    expected = [9.0, 155.0 / 7.0, 214.0 / 7.0, 127.0 / 3.0]
    loo = E.leave_one_out_correction(D, M, line, E.CORRECT_WEIGHT_NONE)
    for index in range(4):
        t.close(loo["predictions"][index], expected[index], "line leave-one-out: field " + str(index + 1))
    t.close(loo["rmse"], 2.576494, "line leave-one-out RMSE", 1e-6)
    t.close(loo["mae"], (3.0 + 22.0 / 7 + 17.0 / 7 + 4.0 / 3) / 4.0, "line leave-one-out MAE")
    # Raw counts: D - M = -2 1 -3 -1, RMSE = sqrt(15/4) = 1.936492, smaller than
    # 2.576 (line) and than the proportional leave-one-out RMSE 1.976932:
    # proportional folds (sum D M = 3640 - D M) / (3000 - D D) x D give
    # 10.379310 21.153846 30.571429 42.571429.
    # Bland-Altman: mean -1.25, SD sqrt(8.75/3) = 1.707825, limits -4.597337 and 2.097337.
    calibration = E.build_count_calibration(E.correct_pairs_from_counts(D, M), line)
    t.close(calibration["agreement"]["rmse"], math.sqrt(15.0 / 4.0), "raw RMSE")
    t.close(calibration["loo"][origin]["rmse"], 1.976932, "proportional leave-one-out RMSE", 1e-6)
    t.close(calibration["agreement"]["loa_low"], -4.597337, "lower limit of agreement", 1e-6)
    t.close(calibration["agreement"]["loa_high"], 2.097337, "upper limit of agreement", 1e-6)
    t.equal((calibration["recommended"][line], calibration["recommended"][origin]), (False, False),
            "neither model beats the raw counts: no correction is recommended")
    t.check(calibration["warnings"][0].startswith("Correction NOT recommended"), "the selected line model's verdict comes first")
    rows = E.correct_fit_rows(calibration, "run")
    t.equal([row["Correction_Recommended"] for row in rows], [0, 0], "the fit table says so for both models")
    t.equal([row["Selected"] for row in rows], [0, 1], "the line model is marked as the applied one")
    t.close(rows[1]["LeaveOneOut_RMSE_Corrected"], 2.576494, "the fit table carries the line leave-one-out RMSE", 1e-6)
    t.equal(rows[1]["Intercept_CI_Includes_Zero"], 1, "the intercept interval holds zero")


def test_correct_zero_objects(E, t, tmp):
    """Fields and ROIs without objects: the proportional model accepts empty
    calibration fields, leave-one-out becomes undefined when the only field
    with objects is left out, and empty ROIs of a real analysis table get a
    corrected count of zero with an upper limit."""
    origin = E.CORRECT_MODEL_ORIGIN
    line = E.CORRECT_MODEL_LINE
    # Calibration D = 0 0 10, M = 0 1 12 (two empty fields, one with a missed object).
    # Proportional: slope = sum D M / sum D D = 120 / 100 = 1.2;
    #   residuals 0, 1, 0 -> SSE 1 on 2 df, SD = sqrt(0.5) = 0.707107
    # Line: mean D 10/3, mean M 13/3, Sxx 600/9, Sxy 690/9 -> slope 1.15,
    #   intercept 13/3 - 11.5/3 = 0.5; residuals -0.5 0.5 0 -> SSE 0.5 on 1 df
    # Leave-one-out: without field 3 every remaining field has 0 detected, so
    # neither model can be fitted and the leave-one-out error is undefined.
    pairs = E.correct_pairs_from_counts([0, 0, 10], [0, 1, 12])
    calibration = E.build_count_calibration(pairs, channel="1", channel_name="DAPI")
    fit = calibration["fit"]
    t.close(fit["slope"], 1.2, "empty fields: proportional slope")
    t.close(fit["residual_sd"], math.sqrt(0.5), "empty fields: residual SD")
    t.equal((fit["d_min"], fit["d_max"]), (0.0, 10.0), "empty fields: the calibrated range starts at zero")
    t.close(calibration["fits"][line]["slope"], 1.15, "empty fields: line slope")
    t.close(calibration["fits"][line]["intercept"], 0.5, "empty fields: line intercept")
    t.close(calibration["fits"][line]["residual_sd"], math.sqrt(0.5), "empty fields: line residual SD")
    loo = calibration["loo"][origin]
    t.equal(loo["predictions"][2], None, "leaving out the only field with objects leaves nothing to fit")
    t.close(loo["predictions"][0], 0.0, "leave-one-out of an empty field: 1.2 x 0 = 0")
    t.close(loo["predictions"][1], 0.0, "leave-one-out of the second empty field")
    t.check(loo["rmse"] != loo["rmse"], "the leave-one-out RMSE is undefined")
    t.equal(calibration["recommended"][origin], None, "no recommendation either way")
    t.check("could not be computed" in calibration["recommendation_texts"][origin], "and that is said")
    t.equal(len([text for text in calibration["warnings"] if text.startswith("Correction NOT recommended")]), 0,
            "an undefined comparison is not reported as a failed one")
    # Bland-Altman of D - M = 0 -1 -2: mean -1, SD 1, limits -2.96 and 0.96;
    # MAPE only over the fields with manual counts: (100 + 2/12 x 100) / 2 = 58.333333
    agreement = calibration["agreement"]
    t.close(agreement["loa_low"], -2.96, "empty fields: lower limit of agreement")
    t.close(agreement["loa_high"], 0.96, "empty fields: upper limit of agreement")
    t.close(agreement["mape"], 175.0 / 3.0, "empty fields: MAPE leaves out the field with manual 0")
    t.equal(agreement["mape_fields"], 2, "two fields enter the MAPE")
    t.check(correct_test_raises(lambda: E.build_count_calibration(pairs, origin, E.CORRECT_WEIGHT_COUNT)),
            "weights 1 / detected refuse a calibration with empty fields")

    # The calibration tables leave the undefined numbers empty.
    out = E.unique_dir(tmp, "correct_zero")
    written = E.write_count_calibration(calibration, out)
    pair_table = E.correct_read_csv(written["pairs_path"])
    t.equal(correct_test_column(pair_table, "LeaveOneOut_Corrected_Count")[2], "", "no leave-one-out value for field 3")
    percent = correct_test_column(pair_table, "Percent_Error")
    t.equal(percent[0], "", "no percent error for a manual count of 0")
    t.close(float(percent[1]), -100.0, "percent error (D - M) / M x 100 of field 2")
    t.close(float(percent[2]), -200.0 / 12.0, "percent error of field 3")
    fit_table = E.correct_read_csv(written["fit_path"])
    t.equal(correct_test_column(fit_table, "Correction_Recommended"), ["", ""], "no recommendation in the fit table")
    t.equal(correct_test_column(fit_table, "LeaveOneOut_RMSE_Corrected"), ["", ""], "undefined RMSE written as empty")

    # A real analysis with two ROIs that hold no object: an oval ROI and a
    # corner triangle that touches two image borders. The image has two disks
    # far from both ROIs, at (200,180) and (220,150), so every count is 0 by
    # construction. (An image without any bright pixel is not used: with a
    # Manual threshold above the image maximum the engine refuses the image;
    # reported to the engine owner.)
    imp = make_image("no objects.tif", slices=1, disks_per_frame=[[(200, 180), (220, 150)]])
    rois = [correct_test_named(OvalRoi(20, 20, 100, 80), "oval"), correct_test_corner()]
    try:
        result = run_image(E, imp, [nuclei_config(E), marker_config(E)], rois, os.path.join(tmp, "correct_zero_runs"))
    finally:
        imp.close()
    t.equal([row["C1_DAPI_Object_Count"] for row in result["summary_rows"]], [0, 0], "an image without disks has no objects")
    done = E.correct_count_table(os.path.join(result["run_dir"], "roi_summary.csv"), calibration, E.unique_dir(tmp, "correct_zero_out"))
    names = done["columns"]
    # 0 detected: estimate 1.2 x 0 = 0; the half width is t2 x SD x sqrt(1 + 0)
    #   = 4.302653 x 0.707107 = 3.042435. 0 is inside the range 0 to 10.
    for row in done["rows"]:
        t.equal((row[names["count"]], row[names["count_low"]]), (0.0, 0.0), "an empty ROI is corrected to 0 (" + row["Parent_ROI_Name"] + ")")
        t.close(row[names["count_high"]], 3.042435, "an empty ROI keeps an upper limit (" + row["Parent_ROI_Name"] + ")", 1e-6)
        t.close(row[names["density_high"]], 3.042435 / float(row["ROI_Area"]), "its density upper limit (" + row["Parent_ROI_Name"] + ")", 1e-6)
        t.equal((row[names["extrapolated"]], row[names["clipped"]]), (0, 0), "0 detected is calibrated here and nothing is clipped")
    t.equal((done["corrected"], done["extrapolated"], done["clipped"], done["blank"]), (2, 0, 0, 0), "row totals for the empty image")


def test_correct_calibrated_frames_and_batch(E, t, tmp):
    """A table from an image calibrated at 0.25 micron per pixel (unit
    'micron') with two frames and three parent ROIs, two of them on the image
    border, and a batch_roi_summary.csv of two images; plus the profile round
    trip and byte-identical analysis tables with the settings changed."""
    root = E.unique_dir(tmp, "correct_frames")
    frame_one = [(30, 30), (80, 30), (200, 30), (30, 150), (200, 150)]
    frame_two = [(30, 30), (200, 30), (200, 150)]

    def calibrated(title, disks_per_frame):
        imp = make_image(title, slices=1, frames=len(disks_per_frame), disks_per_frame=disks_per_frame)
        cal = imp.getCalibration()
        cal.pixelWidth = 0.25
        cal.pixelHeight = 0.25
        cal.setUnit("micron")
        return imp

    def parent_rois():
        # corner: the triangle (0,0) (140,0) (0,140) on the top and left borders.
        #   A disk of radius 6 is inside when x + y < 140 - 6 sqrt(2) = 131.5:
        #   (30,30) and (80,30) are; (200,30), (30,150), (200,150) are far outside.
        # oval: x 150..240, y 100..200 (touches the right and bottom borders),
        #   centre (195,150): holds only the disk at (200,150).
        return [correct_test_named(Roi(0, 0, WIDTH, HEIGHT), "whole"), correct_test_corner(),
                correct_test_named(OvalRoi(150, 100, 90, 100), "oval")]

    # Detected objects by construction, per (frame, ROI):
    expected = {(1, "whole"): 5, (1, "corner"): 2, (1, "oval"): 1, (2, "whole"): 3, (2, "corner"): 1, (2, "oval"): 1}
    configs = [nuclei_config(E), marker_config(E)]
    imp = calibrated("frames micron.tif", [frame_one, frame_two])
    try:
        plain = run_image(E, imp, configs, parent_rois(), os.path.join(root, "plain"), frames=[1, 2])
        changed = run_image(E, imp, configs, parent_rois(), os.path.join(root, "changed"), frames=[1, 2],
                            options_override={"correct_model": E.CORRECT_MODEL_LINE, "correct_weighting": E.CORRECT_WEIGHT_COUNT})
    finally:
        imp.close()
    counts = dict([((row["Frame"], row["Parent_ROI_Name"]), row["C1_DAPI_Object_Count"]) for row in plain["summary_rows"]])
    t.equal(counts, expected, "the two-frame image holds the counts built into it")
    whole = [row for row in plain["summary_rows"] if row["Parent_ROI_Name"] == "whole"][0]
    t.close(whole["ROI_Area"], 240 * 200 * 0.25 * 0.25, "whole-image area is 3000 square microns")

    # With the correction settings changed, the analysis tables are byte for
    # byte the same (apart from the Run_ID, which names the run folder).
    for name in ["roi_summary.csv", "objects.csv", "thresholds.csv"]:
        first = read_utf8(os.path.join(plain["run_dir"], name)).replace(os.path.basename(plain["run_dir"]), "RUN")
        second = read_utf8(os.path.join(changed["run_dir"], name)).replace(os.path.basename(changed["run_dir"]), "RUN")
        t.equal(first, second, name + " does not depend on the correction settings")

    # Correct it with the worked example of this file: slope 6190 / 5500,
    # SD = 1.964457 on 4 df, range 10 to 50. Every count here (1 to 5) is below
    # the range, so every row is flagged as extrapolated.
    calibration = E.build_count_calibration(E.correct_pairs_from_counts(CORRECT_TEST_D, CORRECT_TEST_M), channel="1", channel_name="DAPI")
    summary_path = os.path.join(plain["run_dir"], "roi_summary.csv")
    done = E.correct_count_table(summary_path, calibration, E.unique_dir(root, "out"))
    names = done["columns"]
    slope = 6190.0 / 5500.0
    sd0 = math.sqrt((6982.0 - 6190.0 * 6190.0 / 5500.0) / 4.0)
    table = E.correct_read_csv(done["path"])
    original = E.correct_read_csv(summary_path)
    t.equal(table["headers"][:len(original["headers"])], original["headers"], "original columns first, in order")
    t.equal([row[:len(original["headers"])] for row in table["rows"]], original["rows"], "every original value unchanged")
    for row in done["rows"]:
        key = (int(row["Frame"]), row["Parent_ROI_Name"])
        d = float(expected[key])
        label = " (frame " + str(key[0]) + ", " + key[1] + ")"
        # e.g. 5 detected: 5.627273, half width t4 x SD x sqrt(1 + 25/5500) = 5.466589
        t.close(row[names["count"]], slope * d, "corrected count" + label)
        half = CORRECT_TEST_T4 * sd0 * math.sqrt(1.0 + d * d / 5500.0)
        t.close(row[names["count_high"]], slope * d + half, "upper limit" + label)
        # 5 detected: 5.627273 - 5.466589 = 0.160684; 3 or fewer: below zero, shown as 0
        t.close(row[names["count_low"]], max(0.0, slope * d - half), "lower limit, not below zero" + label)
        t.close(row[names["density"]], slope * d / float(row["ROI_Area"]), "density per square micron" + label)
        t.equal(row[names["extrapolated"]], 1, "below the calibrated range 10 to 50" + label)
    t.close(done["rows"][[r["Parent_ROI_Name"] for r in done["rows"]].index("whole")][names["count_high"]], 5.627273 + 5.466589,
            "upper limit of the whole image in frame 1 as a number", 1e-6)
    t.equal(done["extrapolated"], 6, "all six rows are extrapolated")
    t.equal(len([text for text in done["warnings"] if "more objects than any" in text]), 0,
            "counts below the range do not get the note about larger ROIs")

    # --- a batch of two single-frame images, corrected with real validation runs ---
    folder = os.path.join(root, "batch_in")
    os.makedirs(folder)
    for title, disks in [("p.tif", frame_one), ("q.tif", frame_two)]:
        image = calibrated(title, [disks])
        FileSaver(image).saveAsTiff(os.path.join(folder, title))
        image.close()
    options = E.default_global_options()
    options["save_prepared"] = False
    options["save_overlays"] = False
    batch = E.run_batch(folder, root, configs, configs, options,
                        {"extensions": "tif", "recursive": False, "roi_policy": E.BATCH_ROI_WHOLE})
    t.equal(batch["analyzed"], 2, "both batch images are analyzed")
    batch_table = os.path.join(batch["batch_dir"], "batch_roi_summary.csv")
    cal_dir = os.path.join(root, "calibration")
    os.makedirs(cal_dir)
    correct_test_calibration_runs(E, cal_dir)
    # Pairs (4,6) (5,6) (4,7) (6,9) (7,8) (6,10): slope 252/178, SD 1.359114 on 5 df,
    # range 4 to 7 (see the test with validation runs above).
    # p.tif: 5 detected -> 7.078652, half width t5 x SD x sqrt(1 + 25/178) = 3.731001, inside;
    # q.tif: 3 detected -> 4.247191, half width 3.580949, below the range.
    result = E.run_count_correction(calibration_folder=cal_dir, table_path=batch_table, output_root=os.path.join(root, "batch_out"))
    rows = result["corrected"]["rows"]
    names = result["corrected"]["columns"]
    t.equal([os.path.basename(row["Source_Image"]) if "Source_Image" in row else "" for row in rows], ["p.tif", "q.tif"],
            "one corrected row per batch image")
    t.close(rows[0][names["count"]], 7.078652, "batch: corrected count of p.tif", 1e-6)
    t.close(rows[0][names["count_high"]] - rows[0][names["count"]], 3.731001, "batch: half width for p.tif", 1e-6)
    t.close(rows[1][names["count"]], 4.247191, "batch: corrected count of q.tif", 1e-6)
    t.close(rows[1][names["count"]] - rows[1][names["count_low"]], 3.580949, "batch: half width for q.tif", 1e-6)
    t.equal([row[names["extrapolated"]] for row in rows], [0, 1], "5 is inside the range 4 to 7, 3 is not")
    t.equal(result["corrected"]["settings_match"], "SAME", "the batch settings.properties matches the validation runs")
    t.check(os.path.basename(result["corrected_path"]) == "batch_roi_summary_corrected.csv", "the corrected batch table is named after it")
    t.check(not os.path.exists(os.path.join(batch["batch_dir"], "batch_roi_summary_corrected.csv")), "nothing is written into the batch folder")

    # --- a table path that does not exist creates nothing ---
    missing_root = os.path.join(root, "missing_out")
    t.check(correct_test_raises(lambda: E.run_count_correction(calibration_folder=cal_dir, table_path=os.path.join(root, "nowhere", "roi_summary.csv"),
                                                               output_root=missing_root)), "a missing table is refused")
    t.check(not os.path.exists(missing_root) and not os.path.exists(os.path.join(root, "nowhere")), "and no folder is made for it")

    # --- a very long output path is pointed out (Windows programs cannot open it) ---
    deep = os.path.join(root, "d" * 120, "e" * 120)
    long_result = E.run_count_correction(calibration_folder=cal_dir, table_path=batch_table, output_root=deep)
    t.check(len(long_result["corrected_path"]) > 259, "the test path is longer than 259 characters")
    t.check(len([text for text in long_result["warnings"] if "characters long" in text]) == 1, "and the long path is warned about")
    t.equal(len([text for text in result["warnings"] if "characters long" in text]), 0, "a short path is not")

    # --- profile round trip through write_settings_properties / load_profile_settings ---
    profile = os.path.join(root, "correct_round_trip.properties")
    saved = E.default_global_options()
    saved["correct_model"] = E.CORRECT_MODEL_LINE
    saved["correct_weighting"] = E.CORRECT_WEIGHT_COUNT
    E.write_settings_properties(profile, saved, configs)
    loaded = E.load_profile_settings(profile)[0]
    t.equal((loaded["correct_model"], loaded["correct_weighting"]), (E.CORRECT_MODEL_LINE, E.CORRECT_WEIGHT_COUNT),
            "a profile written with the correction settings loads back with the same values")


def test_correct_large_and_awkward_csv(E, t, tmp):
    """A batch table of realistic size is read quickly and exactly, and line
    breaks inside quoted fields (Unix and Windows) survive."""
    import time
    folder = E.unique_dir(tmp, "correct_big_csv")
    # 2000 rows x 80 columns written by write_csv: about 2 MB. Before the
    # line-by-line parser, a 0.75 MB table took 47 s under Jython and a
    # 5.7 MB batch table did not finish in 15 minutes.
    headers = ["Source_Image", "ROI_Area"] + ["Col_%02d" % index for index in range(77)] + ["C1_DAPI_Object_Count"]
    rows = []
    for index in range(2000):
        row = {"Source_Image": "img " + MICRO + str(index), "ROI_Area": 1000.0 + index, "C1_DAPI_Object_Count": index % 60}
        for column in range(77):
            row["Col_%02d" % column] = 0.123456789 * (index + column)
        rows.append(row)
    path = os.path.join(folder, "batch_roi_summary.csv")
    E.write_csv(path, headers, rows)
    started = time.time()
    table = E.correct_read_csv(path)
    seconds = time.time() - started
    t.check(seconds < 60.0, "a 2 MB table is read in well under a minute (took %.1f s)" % seconds)
    t.equal(len(table["rows"]), 2000, "every row of the large table is read")
    t.equal(table["rows"][1234][0], "img " + MICRO + "1234", "text of a late row")
    t.equal(table["rows"][1999][-1], "19", "last field of the last row (1999 mod 60)")
    t.close(float(table["rows"][777][5]), 0.123456789 * (777 + 3), "a number of a middle row (column Col_03)")
    calibration = E.build_count_calibration(E.correct_pairs_from_counts(CORRECT_TEST_D, CORRECT_TEST_M))
    started = time.time()
    done = E.correct_count_table(path, calibration, E.unique_dir(folder, "out"))
    seconds = time.time() - started
    t.check(seconds < 120.0, "the large table is corrected in reasonable time (took %.1f s)" % seconds)
    t.equal(done["corrected"], 2000, "every row is corrected")

    # Awkward but legal CSV: quoted line breaks (LF and CRLF), a quoted field
    # over three lines, doubled quotes, a literal quote inside an unquoted
    # field, Windows line ends between records, and an unclosed quote at the end.
    text = (u'A,B,C\r\n"one\r\ntwo",x,"say ""hi"""\r\n5" wide,"l1\nl2\nl3",z\r\n'
            u'last,"never closed\nstill open')
    records = E.correct_parse_csv(text)
    t.equal(records[0], ["A", "B", "C"], "header with a Windows line end")
    t.equal(records[1], ["one\r\ntwo", "x", 'say "hi"'], "a CRLF inside quotes is kept; doubled quotes are one quote")
    t.equal(records[2], ['5" wide', "l1\nl2\nl3", "z"], "a quote inside an unquoted field is literal; a field over three lines")
    t.equal(records[3], ["last", "never closed\nstill open"], "an unclosed quote runs to the end of the text")
    t.equal(len(records), 4, "four records")
    t.equal(E.correct_parse_csv(u"a;b\r1;2\r"), [["a", "b"], ["1", "2"]], "old Mac line ends")
    t.equal(E.correct_parse_csv(u"a,b,\n1,,\n"), [["a", "b", ""], ["1", "", ""]], "empty trailing fields are kept")


TESTS.append(("corrected counts: exact data and the t quantile", test_correct_exact_fits))
TESTS.append(("corrected counts: example worked by hand", test_correct_hand_worked_example))
TESTS.append(("corrected counts: refusals, flags and warnings", test_correct_refusals_and_flags))
TESTS.append(("corrected counts: reading tables", test_correct_reads_csv))
TESTS.append(("corrected counts: from validation runs to a corrected table", test_correct_from_validation_runs))
TESTS.append(("corrected counts: settings checks", test_correct_settings_checks))
TESTS.append(("corrected counts: nothing else changes", test_correct_changes_nothing_else))
TESTS.append(("corrected counts: weighted line and line leave-one-out", test_correct_weighted_line_and_line_leave_one_out))
TESTS.append(("corrected counts: fields and ROIs without objects", test_correct_zero_objects))
TESTS.append(("corrected counts: calibrated two-frame table and a batch", test_correct_calibrated_frames_and_batch))
TESTS.append(("corrected counts: large and awkward CSV tables", test_correct_large_and_awkward_csv))
