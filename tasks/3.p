#               **    DATA READING    **

#filename to be read
filename='3'

#Read multiple
#start="3_"
#filelist=system("ls ".start."*.csv")
#print filelist

#For using ',' as a decimal sign
#set decimalsign locale

#For csv files
#set datafile separator ","

#               **      AXIS SETTINGS AND TITLE        **

#axis labels
set ylabel "Efficiency"
set xlabel "Energy (keV)"

#set logscale x
#set logscale y

#graph title
#set title "Diffrakció I-P karakterisztikája"

#Scientifix notation the way I like it
#set format y "%2.0t{/Symbol \264}10^{%L}"

#               **      OUTPUT FORMAT   **

#eps
set terminal eps
set output filename.'.eps'

#png
#set terminal png size 960,540
#set output filename.'.png'

#              **      FUNCTION FITTING        **

#f(x)=b+sin(a*sqrt(x))**2
#a = 1
#b=1
#fit[0:0.4] f(x) filename u 1:2:3:4 xyerr via a,b

#               **      PLOTTING        **

#crop a subplot to a smaller range
#u ($1>1.7e-7?$1:1/0):($2)

#skip first 5 lines
#every ::5

#set key autotitle columnhead

#change plot titles size
#set key font "0,5"

#error box
#boxxyerr

#plot single
p filename using 1:2 t "η_t_o_t" lw 2, filename using 1:4 t "η_i_n_t" lw 2

#plot multiple
#p for [fn in filelist] fn u ($1):(sqrt(10**(($2)/10))) w lines lw 2
