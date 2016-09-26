rmdir "../ProbeStructures/test" /s /q
mkdir "../ProbeStructures/test"
xcopy "../ProbeStructures/currenter" "../ProbeStructures/test" /e

rmdir "../tetgen/test" /s /q
mkdir "../tetgen/test"
xcopy "../tetgen/current" "../tetgen/test" /e