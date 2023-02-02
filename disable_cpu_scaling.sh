cpu_num=$(cat /proc/cpuinfo | grep "processor" | grep ":" | wc -l)

for ((ix = 0; ix < ${cpu_num}; ix++)); do
    echo "performance" >/sys/devices/system/cpu/cpu${ix}/cpufreq/scaling_governor
done