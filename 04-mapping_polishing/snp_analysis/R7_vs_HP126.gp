set terminal png tiny size 800,800
set output "R7_vs_HP126.png"
set xtics rotate ( \
 "contig_2_pilon" 1.0, \
 "contig_3_pilon" 9361213.0, \
 "" 9653782 \
)
set ytics ( \
 "contig_3_pilon" 1.0, \
 "contig_4_pilon" 551066.0, \
 "contig_1_pilon" 1605450.0, \
 "*contig_2_pilon" 8601362.0, \
 "" 8771744 \
)
set size 1,1
set grid
unset key
set border 0
set tics scale 0
set xlabel "REF"
set ylabel "QRY"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
set xrange [1:9653782]
set yrange [1:8771744]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "R7_vs_HP126.fplot" title "FWD" w lp ls 1, \
 "R7_vs_HP126.rplot" title "REV" w lp ls 2
