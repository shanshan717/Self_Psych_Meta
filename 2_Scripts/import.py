import os

base_dir = "/Volumes/ss/Self_Psych_Meta/3_Output/1_ALE"

for fname in os.listdir(base_dir):
    old_path = os.path.join(base_dir, fname)
    new_name = fname

    # Rule 1: control_minus_patient -> hypoactivation
    if fname.startswith("control_minus_patient"):
        new_name = fname.replace(
            "control_minus_patient",
            "hypoactivation",
            1
        )

    # Rule 2: patient_minus_control -> hyperactivation
    if fname.startswith("patient_minus_control"):
        new_name = fname.replace(
            "patient_minus_control",
            "hyperactivation",
            1
        )

    new_path = os.path.join(base_dir, new_name)

    if new_name != fname:
        os.rename(old_path, new_path)
        print(f"Renamed: {fname} -> {new_name}")