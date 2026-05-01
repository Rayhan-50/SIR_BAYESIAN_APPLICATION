"""
SIR Model Data Preprocessor - Vectorized Version
Converts 4 monthly Excel files to weekly frequency (2025+),
then merges with 2 already-weekly files into unified_weekly_dataset.xlsx
"""

import sys, warnings
import pandas as pd
import numpy as np
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

DATA_DIR = Path(r"d:\research\research\sir_model_application\data")

MONTH_MAP = {
    "jan":1,"feb":2,"mar":3,"apr":4,"may":5,"jun":6,
    "jul":7,"aug":8,"sep":9,"oct":10,"nov":11,"dec":12
}

def banner(t): print(f"\n{'='*65}\n  {t}\n{'='*65}")
def sub(m):    print(f"  >> {m}")

# ── core helpers ──────────────────────────────────────────────────────────────

def week_monday(s):
    """Snap a date Series to the Monday of its ISO week (vectorized)."""
    s = pd.to_datetime(s)
    return s - pd.to_timedelta(s.dt.weekday, unit="D")

def parse_year_month(df, yc, mc):
    """Build a month-start Timestamp from Year (numeric) + Month (text/numeric)."""
    year = pd.to_numeric(df[yc], errors="coerce")
    month = df[mc].apply(lambda m: (
        MONTH_MAP.get(str(m).strip().lower()[:3])
        if isinstance(m, str)
        else (int(float(m)) if pd.notna(m) else None)
    ))
    valid = year.notna() & month.notna() & year.between(2000, 2030) & month.between(1, 12)
    dates = pd.Series(pd.NaT, index=df.index)
    dates[valid] = pd.to_datetime(
        year[valid].astype(int).astype(str) + "-" +
        month[valid].astype(int).astype(str).str.zfill(2) + "-01"
    )
    return dates

def monthly_to_weekly_vectorized(df, date_col, value_cols):
    """
    Vectorized monthly->weekly proportional expansion.
    Repeats each month-row once per ISO week that overlaps that month,
    dividing values by the number of overlapping weeks.
    """
    df = df.copy()
    df[date_col] = pd.to_datetime(df[date_col])
    df = df.dropna(subset=[date_col])

    records = []
    for ms, grp in df.groupby(date_col):
        me = ms + pd.offsets.MonthEnd(0)
        # All days in the month, snap to Monday
        days = pd.date_range(ms, me, freq="D")
        mondays = pd.Series(days - pd.to_timedelta(days.weekday, unit="D")).unique()
        mondays = pd.DatetimeIndex(sorted(mondays))
        n = len(mondays)
        # Sum all rows in this month-group first, then divide
        sums = grp[value_cols].sum(numeric_only=True)
        for mon in mondays:
            row = {"week_start": mon}
            for vc in value_cols:
                row[vc] = sums[vc] / n if vc in sums else np.nan
            records.append(row)

    if not records:
        return pd.DataFrame(columns=["week_start"] + value_cols)

    out = pd.DataFrame(records)
    out["week_start"] = pd.to_datetime(out["week_start"])
    return out.groupby("week_start", as_index=False)[value_cols].sum()

def find_header_row(path, sheet=0):
    """Find the row index that contains 'Year' or 'year'."""
    raw = pd.read_excel(path, sheet_name=sheet, header=None, nrows=15)
    for i, row in raw.iterrows():
        if any(str(v).strip().lower() in ("year","yr") for v in row if pd.notna(v)):
            return i
    return 0

def load_nisb(path, label):
    """Generic loader for NISB-style files with Year + Month columns."""
    banner(label)
    hr = find_header_row(path)
    df = pd.read_excel(path, sheet_name=0, header=hr)
    sub(f"Shape: {df.shape}  |  Cols: {list(df.columns)}")

    yc = next((c for c in df.columns if str(c).strip().lower() in ("year","yr")), None)
    mc = next((c for c in df.columns if str(c).strip().lower() in ("month","mo")), None)
    sub(f"Year col: {yc!r}   Month col: {mc!r}")

    if not yc or not mc:
        raise ValueError(f"Cannot find Year/Month columns in {label}")

    # Forward-fill merged year cells
    df[yc] = pd.to_numeric(df[yc], errors="coerce").ffill()
    df["_date"] = parse_year_month(df, yc, mc)
    df = df.dropna(subset=["_date"]).copy()

    skip = {yc, mc, "_date"}
    vcols = [c for c in df.columns
             if c not in skip
             and pd.api.types.is_numeric_dtype(df[c])
             and df[c].notna().any()]
    sub(f"Value cols ({len(vcols)}): {vcols}")

    df_all_years = df.copy()
    df_2025 = df[df["_date"].dt.year >= 2025].copy()
    sub(f"Rows >= 2025: {len(df_2025)}  (total: {len(df)})")
    if df_2025.empty:
        sub("WARNING: No 2025+ rows — using ALL years as fallback")
        df_2025 = df_all_years.copy()

    weekly = monthly_to_weekly_vectorized(df_2025, "_date", vcols)
    sub(f"Weekly rows: {len(weekly)}")
    if not weekly.empty:
        sub(f"Date range: {weekly['week_start'].min().date()} -> {weekly['week_start'].max().date()}")
    return weekly

def load_holidays(path):
    banner("FILE 4: Govt School Holidays")
    hr = find_header_row(path)
    df = pd.read_excel(path, sheet_name=0, header=hr)
    sub(f"Shape: {df.shape}  |  Cols: {list(df.columns)}")
    sub(f"First 5 rows:\n{df.head().to_string()}\n")

    # Find a proper Date column (exact daily dates, not Year which looks date-like)
    date_col = None
    for col in df.columns:
        if str(col).strip().lower() in ("year", "yr"):
            continue  # skip Year column
        parsed = pd.to_datetime(df[col], errors="coerce")
        if parsed.notna().mean() > 0.5:
            df[col] = parsed
            date_col = col
            break
    sub(f"Daily date column: {date_col!r}")

    # Find holiday value column (0/1 or counts)
    val_col = next((c for c in df.columns
                    if str(c).strip().lower() in ("holiday", "is_holiday", "flag", "count", "days")
                    and c != date_col), None)
    if val_col is None:
        # fallback: first numeric col that is not Year
        for c in df.columns:
            if c == date_col:
                continue
            if pd.api.types.is_numeric_dtype(df[c]) and str(c).strip().lower() not in ("year","yr"):
                val_col = c
                break
    sub(f"Holiday value column: {val_col!r}")

    if date_col is not None:
        df[date_col] = pd.to_datetime(df[date_col], errors="coerce")
        df_2025 = df[df[date_col].dt.year >= 2025].copy()
        sub(f"Rows >= 2025: {len(df_2025)}")
        if df_2025.empty:
            sub("WARNING: No 2025+ rows — using ALL years")
            df_2025 = df.copy()

        # Vectorised: snap each date to its Monday
        df_2025 = df_2025.dropna(subset=[date_col]).copy()
        df_2025["week_start"] = week_monday(df_2025[date_col])

        if val_col:
            df_2025[val_col] = pd.to_numeric(df_2025[val_col], errors="coerce").fillna(0)
            weekly = df_2025.groupby("week_start", as_index=False)[val_col].sum()
            weekly = weekly.rename(columns={val_col: "holiday_days"})
        else:
            # Count rows per week as proxy
            weekly = df_2025.groupby("week_start").size().reset_index(name="holiday_days")
    else:
        sub("WARNING: No date column found — empty output")
        weekly = pd.DataFrame(columns=["week_start", "holiday_days"])

    weekly["week_start"] = pd.to_datetime(weekly["week_start"])
    weekly = weekly.sort_values("week_start").reset_index(drop=True)
    sub(f"Weekly rows: {len(weekly)}")
    if not weekly.empty:
        sub(f"Date range: {weekly['week_start'].min().date()} -> {weekly['week_start'].max().date()}")
    return weekly

def load_already_weekly(path, label, agg="sum"):
    banner(f"ALREADY-WEEKLY: {label}")
    df = pd.read_excel(path, sheet_name=0, header=0)
    sub(f"Shape: {df.shape}  |  Cols: {list(df.columns)}")
    sub(f"First 5 rows:\n{df.head().to_string()}\n")

    dc = None
    # 1) Only try datetime parsing on string/object dtype columns (avoids int cols like Case~200)
    for col in df.columns:
        if str(col).strip().lower() in ("year","yr","week","wk","epiweek"):
            continue
        if not (df[col].dtype == object or str(df[col].dtype).startswith("datetime")):
            continue
        parsed = pd.to_datetime(df[col], errors="coerce")
        if parsed.notna().mean() > 0.4:
            df[col] = parsed
            dc = col
            break

    # 2) Year + ISO Week Number (integer Week column)
    if dc is None:
        yc  = next((c for c in df.columns if str(c).strip().lower() in ("year","yr")), None)
        wkc = next((c for c in df.columns if str(c).strip().lower() in ("week","wk","epiweek","isoweek")), None)
        if yc and wkc:
            sub(f"Parsing ISO week: Year={yc!r} + Week={wkc!r}")
            df[yc]  = pd.to_numeric(df[yc],  errors="coerce")
            df[wkc] = pd.to_numeric(df[wkc], errors="coerce")
            valid = df[yc].notna() & df[wkc].notna()
            df["_date"] = pd.NaT
            df.loc[valid, "_date"] = df[valid].apply(
                lambda r: pd.to_datetime(
                    f"{int(r[yc])}-W{int(r[wkc]):02d}-1", format="%G-W%V-%u"
                ), axis=1
            )
            dc = "_date"

    # 3) Year + Month
    if dc is None:
        yc = next((c for c in df.columns if str(c).strip().lower() in ("year","yr")), None)
        mc = next((c for c in df.columns if str(c).strip().lower() in ("month","mo")), None)
        if yc and mc:
            df[yc] = pd.to_numeric(df[yc], errors="coerce").ffill()
            df["_date"] = parse_year_month(df, yc, mc)
            dc = "_date"

    # 4) Year-only fallback: assign sequential week numbers within each year
    if dc is None:
        yc = next((c for c in df.columns if str(c).strip().lower() in ("year","yr")), None)
        if yc:
            sub(f"Fallback: Year-only col {yc!r}. Using row rank within year as week number.")
            df[yc] = pd.to_numeric(df[yc], errors="coerce")
            df = df.dropna(subset=[yc]).copy()
            df["_wk"] = df.groupby(yc).cumcount() + 1
            valid = df[yc].between(2000, 2030) & df["_wk"].between(1, 53)
            df["_date"] = pd.NaT
            df.loc[valid, "_date"] = df[valid].apply(
                lambda r: pd.to_datetime(
                    f"{int(r[yc])}-W{int(r['_wk']):02d}-1", format="%G-W%V-%u"
                ), axis=1
            )
            dc = "_date"

    sub(f"Date column used: {dc!r}")
    if dc is None:
        sub("WARNING: No date column — skipping")
        return pd.DataFrame(columns=["week_start"])

    df = df.dropna(subset=[dc]).copy()
    df[dc] = pd.to_datetime(df[dc], errors="coerce")
    df_2025 = df[df[dc].dt.year >= 2025].copy()
    sub(f"Rows >= 2025: {len(df_2025)}")
    if df_2025.empty:
        sub("WARNING: No 2025+ data — using all years")
        df_2025 = df.copy()

    df_2025["week_start"] = week_monday(df_2025[dc])
    skip_cols = {"year","yr","week","wk","epiweek","isoweek","_wk"}
    vcols = [c for c in df_2025.columns
             if c not in (dc, "_date", "_wk", "week_start")
             and str(c).strip().lower() not in skip_cols
             and pd.api.types.is_numeric_dtype(df_2025[c])
             and df_2025[c].notna().any()]
    sub(f"Value cols: {vcols}")

    fn = df_2025.groupby("week_start", as_index=False)[vcols].sum if agg == "sum" else \
         df_2025.groupby("week_start", as_index=False)[vcols].mean
    weekly = fn()
    sub(f"Weekly rows: {len(weekly)}  |  {weekly['week_start'].min().date()} -> {weekly['week_start'].max().date()}")
    return weekly

# ── run all ───────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    results = []

    # File 1: Enrolment
    w = load_nisb(
        DATA_DIR / "(NISB) Monthly Enrolment by Different Age Groups 2022 to 2025.xlsx",
        "FILE 1: Monthly Enrolment"
    )
    out1 = DATA_DIR / "(NISB) Monthly Enrolment by Different Age Groups 2022 to 2025_weekly.xlsx"
    w.to_excel(out1, index=False)
    sub(f"Saved: {out1.name}")
    results.append(("Enrolment_weekly", w.rename(columns=lambda c: f"enrol_{c}" if c != "week_start" else c)))

    # File 2: ILI
    w = load_nisb(
        DATA_DIR / "(NISB) Sentinel Site Wise Influenza like illness (ILI) Lab Results 2022 to 2025.xlsx",
        "FILE 2: Sentinel ILI Lab Results"
    )
    out2 = DATA_DIR / "(NISB) Sentinel Site Wise Influenza like illness (ILI) Lab Results 2022 to 2025_weekly.xlsx"
    w.to_excel(out2, index=False)
    sub(f"Saved: {out2.name}")
    results.append(("ILI_weekly", w.rename(columns=lambda c: f"ili_{c}" if c != "week_start" else c)))

    # File 3: Severe
    w = load_nisb(
        DATA_DIR / "(NISB) Severe case.xlsx",
        "FILE 3: Severe Cases"
    )
    out3 = DATA_DIR / "(NISB) Severe case_weekly.xlsx"
    w.to_excel(out3, index=False)
    sub(f"Saved: {out3.name}")
    results.append(("Severe_weekly", w.rename(columns=lambda c: f"severe_{c}" if c != "week_start" else c)))

    # File 4: Holidays
    w = load_holidays(DATA_DIR / "Govt_School_Holiday_No_Year_2022 to 2025 (1).xlsx")
    out4 = DATA_DIR / "Govt_School_Holiday_No_Year_2022 to 2025 (1)_weekly.xlsx"
    w.to_excel(out4, index=False)
    sub(f"Saved: {out4.name}")
    results.append(("Holiday_weekly", w.rename(columns=lambda c: f"holiday_{c}" if c != "week_start" else c)))

    # Already-weekly: pct positive
    w = load_already_weekly(
        DATA_DIR / "% percentage of specimens positive for influenza 1.xlsx",
        "Pct Positive Influenza", agg="mean"
    )
    results.append(("PctPositive_weekly", w.rename(columns=lambda c: f"pct_{c}" if c != "week_start" else c)))

    # Already-weekly: SARI
    w = load_already_weekly(
        DATA_DIR / "Number of SARI weekly cases 1.xlsx",
        "SARI Weekly Cases", agg="sum"
    )
    results.append(("SARI_weekly", w.rename(columns=lambda c: f"sari_{c}" if c != "week_start" else c)))

    # Merge all
    banner("MERGING ALL")
    merged = results[0][1]
    for name, df in results[1:]:
        if df.empty or "week_start" not in df.columns:
            sub(f"Skipping empty: {name}")
            continue
        merged = pd.merge(merged, df, on="week_start", how="outer")
    merged = merged.sort_values("week_start").reset_index(drop=True)
    merged.to_excel(DATA_DIR / "unified_weekly_dataset.xlsx", index=False)
    sub(f"Unified shape: {merged.shape}")
    sub(f"Saved: unified_weekly_dataset.xlsx")

    # Summary table
    banner("SUMMARY TABLE")
    print(f"\n  {'Source':<45} {'Rows':>5}  {'From':<12}  {'To':<12}")
    print(f"  {'-'*45}  {'-'*5}  {'-'*12}  {'-'*12}")
    for name, df in results:
        ws = df.get("week_start") if hasattr(df, "get") else df["week_start"] if "week_start" in df.columns else None
        if ws is not None and ws.notna().any():
            s, e = str(ws.min())[:10], str(ws.max())[:10]
        else:
            s = e = "N/A"
        print(f"  {name:<45} {len(df):>5}  {s:<12}  {e:<12}")
    print(f"  {'UNIFIED':<45} {len(merged):>5}  {str(merged['week_start'].min())[:10]:<12}  {str(merged['week_start'].max())[:10]:<12}")

    banner("ALL DONE")
    print("\n  Saved files:")
    for f in sorted(DATA_DIR.glob("*_weekly.xlsx")):
        print(f"    {f.name}")
    print("    unified_weekly_dataset.xlsx")
