#!/usr/bin/env python3
import sys
import os
import shutil

def split_maf(file_in, dir_out, max_num):
    # Ensure output directory exists
    os.makedirs(dir_out, exist_ok=True)

    # Read all lines
    with open(file_in, 'r') as f:
        lines = f.readlines()

    # Separate header and data lines
    header = lines[0]
    data_lines = lines[1:]
    total_lines = len(data_lines)

    # If fewer data lines than splits, just copy file
    if total_lines < max_num:
        base_name = os.path.basename(file_in)
        out_path = os.path.join(dir_out, f"{base_name}.1")
        shutil.copyfile(file_in, out_path)
        print(f"Copied original file to {out_path} (less than {max_num} lines)")
        return

    # Determine lines per file
    lines_per_file = total_lines // max_num
    for i in range(max_num):
        start = i * lines_per_file
        # last file gets any remaining lines
        end = (i + 1) * lines_per_file if i < max_num - 1 else total_lines

        out_lines = [header] + data_lines[start:end]

        base_name = os.path.basename(file_in)
        out_file = os.path.join(dir_out, f"{base_name}.{i+1}")
        with open(out_file, 'w') as f_out:
            f_out.writelines(out_lines)

        print(f"Wrote {len(out_lines) - 1} lines to {out_file}")

# --------- Script entry point ----------
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python split_maf.py <file_in> <dir_out> <max_num>")
        sys.exit(1)

    file_in = sys.argv[1]
    dir_out = sys.argv[2]
    max_num = int(sys.argv[3])

    split_maf(file_in, dir_out, max_num)

