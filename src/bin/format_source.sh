#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

files=(`find ./src -name '*.h' -o -name '*.hpp' -o -name '*.c' -o -name '*.cpp' -o -name '*.glsl'`)

filecount=${#files[@]}

count=0
for f in ${files[@]}
do
	# Run it in parallel
	echo -ne Formatting ${count}/${filecount}\\r
	${DIR}/clang-format -sort-includes=false -i ${f} &
	count=$((${count}+1))

	# Throttle the parallel commands
	if !((count % 16)); then
		wait
	fi
done

wait

echo Done! Formatted ${filecount} files


