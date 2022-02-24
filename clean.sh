#! /usr/bin/env bash

# delete the files and symlinks (but don't follow links)
if [[ $(which munlink) ]]; then
    find -P work/?? -type f -print0 -o -type l -print0 | xargs -0 munlink 2&>/dev/null
else
    find -P work/?? -type f -delete -o -type l -delete
fi

# delete the empty directories
find -P work/?? -type d -empty -delete

# delete the back up reporting files
rm dag.dot.* report.html.* timeline.html.* trace.txt.*
