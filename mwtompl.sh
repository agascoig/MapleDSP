#!/bin/sh -e

maple=/Library/Frameworks/Maple.framework/Versions/2026/bin/maple

for file in *.mw; do
# Skip if no matches (zsh leaves literal if none unless nullglob is set)
   [[ -e "$file" ]] || continue

   out="${file%.mw}.mpl" # remove .mw extension, add .mpl

   echo "Checking $file -> $out"

# Get modification times (seconds since epoch)
    if [[ -e "$out" ]]; then
        mw_time=$(stat -f %m "$file")
        mpl_time=$(stat -f %m "$out")
        
        # Difference in seconds
        diff=$(( mw_time - mpl_time ))
        
        if [[ $diff -gt 10 ]]; then
            echo "  → Converting (source is newer by ${diff}s)"
            echo "a:=Worksheet:-WorksheetToMapleText(\"${file}\"):printf(\"%s\\\n\", a);" | "$maple" -q > "${out}"
        else
            echo "  → Skipping (up to date)"
        fi
    else
        echo "  → Converting (no output file exists)"
        echo "a:=Worksheet:-WorksheetToMapleText(\"${file}\"):printf(\"%s\\\n\", a);" | "$maple" -q > "${out}"
    fi
done
