Load in the pickle and then run the following commands to generate the new Figure 12b plot for the Review Paper:


-= Main Plot: =-

plt.clf();plt.figure();plt.errorbar(plotTemp,plotMob,yerr=plotMobError, color='k');plt.yscale('log');plt.xlabel('Temperature, K');plt.ylabel('Mobility, cm'+r'$^{2}$ '+'V'+r'$^{-1}$'+r's$^{-1}$');plt.xlim([280,500]);plt.ylim([9E-6,2E0]);plt.savefig('./figure3SI.pdf')


-= Exp Data: =-

plt.plot([280,500],[1.21E-3,1.21E-3],'--r',zorder=-32),plt.plot([280,500],[1E-5,1E-5],'--g',zorder=-32);plt.plot([280,500],[1E-3,1E-3],'--g',zorder=-32);plt.plot([280,500],[1E-4,1E-4],'--b',zorder=-32);plt.plot([280,500],[2E-5,2E-5],'--c',zorder=-32);plt.plot([280,500],[4.5E-4,4.5E-4],'--c',zorder=-32);plt.savefig('./figure3SI.pdf')

OR

plt.plot([280,500],[1.21E-3,1.21E-3],'--',color='0.8',zorder=-32),plt.plot([280,500],[1E-5,1E-5],'-.',color='0.2',zorder=-32);plt.plot([280,500],[1E-3,1E-3],'-.',color='0.2',zorder=-32);plt.plot([280,500],[1E-4,1E-4],'--',color='0.6',zorder=-32);plt.plot([280,500],[2E-5,2E-5],':',color='0.4',zorder=-32);plt.plot([280,500],[4.5E-4,4.5E-4],':',color='0.4',zorder=-32);plt.savefig('./figure3SI.pdf')


-= Order Line: =-

plt.plot([317.85,317.85],[1E-5,2],'--m',zorder=-33,linewidth=3);plt.savefig('./figure3SI.pdf')
