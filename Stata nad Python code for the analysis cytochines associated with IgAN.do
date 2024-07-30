clear

**# Import dataset and prepare cytokine datasets for analyses


import excel "C:\Documenti\Sara Albrandi\Olink_IgAN_clean.xlsx", sheet("Olink + baseline") firstrow
cd "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi"

**# prepare the master dataset
preserve
rename tewak tweak
keep group il8-csf1
foreach var of varlist il8-csf1 {
  rename `var' `=strupper("`var'")'
}
cap drop unwanted
egen unwanted = rowmiss(_all)
drop if unwanted
drop unwanted
save IgANcytochines, replace
export delimited _all using IgANcytochines, replace
restore


label var group "Group"
label define group 0 "Controls" 1 "IgA Nephropathy"
label values group group
label define yesno 0 "No" 1 "Yes"
label var gender "Sex"
label define gender 0 "Male" 1 "Female"
label values gender gender
label var race "Ethnicity"
label define race 0 "White" 1 "Black" 2 "Hispanic" 3 "Asian" 4 "Other"
label values race race
label var is "Immunosuppression"
label values is yesno 
label var raas "RASS inhibitors"
label values raas yesno
label var proteinuria "Proteinuria, gr/day"
label var serum_creat "Serum Creatinine, mg/dL"
label var egfr "eGFR, ml/min/1.73m2"
label var serum_alb "Serum Albumin, g/dL"
label var hematuria "Hematuria, RBC/ul" 
* Oxford MEST-C score
cap drop mestc
gen mestc = m + e + s + t + c
label var mestc "MEST-C score"
label var m "m in MEST-C score"
label var s "s in MEST-C score"
label var e "e in MEST-C score"
label var t "t in MEST-C score"
label var c "c in MEST-C score"
label var iga_elisa "Total IgA, ug/ml"
label var gdiga_elisa "gd-IgA1 , ug/ml"
label var iggaa_company2 "IgG autoantibodies against gd-IgA1, ug/ml"

**# Table 1

summ age gender race is raas proteinuria serum_creat egfr serum_alb hematuria mestc m s e t c ///
iga_elisa gdiga_elisa iggaa_company2

dtable age i.gender i.race i.is i.raas proteinuria serum_creat egfr serum_alb hematuria mestc m s e t c ///
iga_elisa gdiga_elisa iggaa_company2 ///
	,   ///	
	by(group, tests) ///
	define(meansd = mean sd, delimiter(" ± ")) ///
	define(myiqr = p25 p75, delimiter("-")) ///
	define(myrange = min max, delimiter("-")) ///
	factor(gender race is raas, test(fisher)) ///
	continuous(age proteinuria egfr serum_alb , stat(meansd) test(kwallis)) ///
    continuous(proteinuria serum_creat hematuria mestc m s e t c iga_elisa gdiga_elisa iggaa_company2, stat(median myiqr) test(kwallis)) ///   
	column(by(hide) test(p-value)) ///
	title(Table 1. Characteristic of the groups ) ///
	note(Mann-Whitney test for continuous variables (reported as mean ± standard deviation, or median (interquartile range)), Fisher's exact test for categorical variable) ///
	sformat("%s" sd) ///
	nformat("%3.2f" mean sd median p25 p75) ///
	nformat("%3.1f" min max) ///
	sformat("(%s)" myiqr myrange) ///
    nformat("%3.0f" N count fvfrequency) ///
    nformat("%3.1f" fvpercent ) ///
    nformat("%6.3f" kwallis fisher) ///
	export(table1_R1.html, replace)
	collect export Table1_R1.xlsx, replace
	collect export Table1_R1.docx, replace
	

**# make univariate test with Holm correction
#delimit ;
global allvars "IL8 VEGFA CD8A MCP3 GDNF CDCP1 CD244 IL7 OPG LAP_TGFB1 
	UPA IL6 IL17C MCP1 IL17A CXCL11 AXIN1 TRAIL IL20RA CXCL9 CST5 IL2RB OSM 
	CXCL1 CCL4 CD6 SCF IL18 SLAMF1 TGFA MCP4 CCL11 TNFSF14 FGF23 IL10RA FGF5 
	MMP1 LIFR FGF21 CCL19 IL15RA IL10RB IL18R1 PDL1 CXCL5 TRANCE HGF IL12B 
	MMP10 IL10 TNF CCL23 CD5 CCL3 FIT3L CXCL6 CXCL10 EBP1 SIRT2 CCL28 DNER 
	ENRAGE CD40 IFNGAMMA FGF19 MCP2 CASP8 CCL25 CX3CL1 TNFRSF9 NT3 TWEAK 
	CCL20 ST1A1 STAMBP ADA TNFB CSF1"
	;
 #delimit cr	

di wordcount("$allvars")


cd "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi"
use IgANcytochines, clear

local n_vars: word count $allvars
matrix stats = J(`n_vars',2,.)

local irow=0
foreach varname in $allvars   {
	local ++irow
	qui kwallis `varname', by(group) 
	matrix stats[`irow',1] = r(chi2_adj)
	matrix stats[`irow',2] = chiprob(r(df), r(chi2_adj))
  }
matrix rownames stats= $allvars
matrix colnames stats= "Chi2" "Pval"
matrix list stats


preserve
clear
 mata:
 rownames = st_matrixrowstripe("stats")[.,2]
 rownames
 st_addobs(rows(rownames))
 idx = st_addvar("str33",("Cytokine"))
 st_sstore(1::rows(rownames),idx,rownames)
 idx = st_addvar("double",("Chi2","Pval"))
 st_store(1::rows(rownames), idx, st_matrix("stats"))
 end
qqvalue Pval, method(holm) qvalue(AdjPval)
gen Sig = cond(AdjPval<0.05,1,0)
label define Sig 0 "No" 1 "Yes"
label values Sig Sig
format Chi2 %3.1f
format Pval %4.3f
format AdjPval %4.3f
list, sepby(Cytokine) table noobs
keep Cytokine Chi2 Pval AdjPval Sig
cap save cytokine_adjusted_holm, replace
cap export delimited _all using cytokine_adjusted_holm, replace
restore

**# Table 2
cd "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi"
use IgANcytochines, clear

dtable $allvars ///
	, novarlabel  ///	
	by(group, tests) ///
	define(meansd = mean sd, delimiter(" ± ")) ///
	define(myiqr = p25 p75, delimiter("-")) ///
	define(myrange = min max, delimiter("-")) ///
    continuous($allvars , stat(median myiqr) test(kwallis)) ///   
	column(by(hide) test(p-value)) ///
	title(Table 1. Difference in cytokines between groups ) ///
	note(Mann-Whitney test for continuous variables (reported as mean ± standard deviation).) ///
	sformat("%s" sd) ///
	nformat("%3.1f" mean sd median p25 p75) ///
	nformat("%3.1f" min max) ///
	sformat("(%s)" myiqr myrange) ///
    nformat("%3.0f" N count fvfrequency) ///
    nformat("%3.1f" fvpercent ) ///
    nformat("%6.3f" kwallis fisher) ///
	export(table2_1.html, replace)
collect export Table2_R1.xlsx, replace
collect export Table2_R1.docx, replace



preserve
import excel "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi\Table2_R1.xlsx", sheet("Sheet1") cellrange(A3:E81) firstrow clear
rename N cytokine
sort cytokine
cap save table2_R1, replace
clear
import delimited "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi\cytokine_adjusted_holm.csv"
sort cytokine
cap save holm, replace
clear
use table2_R1
merge 1:1 cytokine using holm
drop _merge
gsort -chi2
rename B Control
rename C IgAN
rename D Total
rename E Pvalue
drop chi2
drop pval
format adjpval %5.4f
rename adjpval Adjpval
rename sig Significant
export excel _all using "Table2_sorted_R1", firstrow(variables) replace
restore



**# Compare Adaptive Lasso, Lasso CV, Elastic Net among significant vars
#delimit ;
global holmvars	"IL6 CCL4 FGF23 IL15RA PDL1 IL12B MMP10 TNF CCL23 
	CD5 CCL3 CD40 CX3CL1 TNFRSF9 CSF1"
	;
#delimit cr

qui lasso logit group $holmvars, selection(adaptive, step(5)) rseed(123)
estimate store adaptlasso
qui lasso logit group $holmvars, rseed(123)
estimates store autolasso
qui elasticnet logit group $holmvars, rseed(123)
estimates store autoenet
lassogof adaptlasso autolasso autoenet

**# Check best number of steps in adaptive Lasso

foreach num of numlist 5(1)15  {
	qui lasso logit group $holmvars, selection(adaptive, step(`num')) rseed(140867)
	estimate store adaptlasso_`num'
	 }
	 
lassogof adaptlasso_*


**# Fit Adaptive Lasso and plot coefpath
cd "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi"
use IgANcytochines, clear
lasso logit group $holmvars, selection(adaptive, step(9)) rseed(7119) nolog
global vars_sel `e(allvars_sel)'
global selected_lambda = e(ID_cv)
lassoselect id =   $selected_lambda
cvplot
global selected_var = e(allvars_sel)
display "$vars_sel"
lassocoef, display(coef) sort(coef)
coefpath, data(coefpath_data, replace) nodraw

preserve
keep group $selected_var
export delimited _all using selected_var, replace
restore


preserve
use coefpath_data, clear       
lookfor $vars_sel

qui summ l1norm if id == $selected_lambda
local xcut = r(mean)     
lookfor $vars_sel

twoway (line var1 l1norm, lwidth(*1.5)) ///
	(line var6 l1norm, lwidth(*1.5)) ///
	(line var13 l1norm, lwidth(*1.5)) ///
	if alpha == 1 & step == 2, yline(0, lcolor("scheme foreground") lpattern(solid)) ///
	ytitle(Standardized coefficients) title("LASSO: Coefficient Path") ylabel(,format(%3.1f)) /// 
	xline(`xcut', lpattern(dash) lcolor(stc2) lwidth(*1.2)) ///
	legend(on)
graph export FigureS1.png, replace
graph export FigureS1.pdf, replace
restore


**# Matrix graph ggally (PairGrid from seaborn) in Python
cd "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi"
python:
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
# change font
plt.rcParams['font.family'] = "Arial"
df = pd.read_csv("selected_vars.csv")
df.dropna(inplace=True)

group_labels = {0: "Controls", 1: "IgA Nephropathy",}
df['group'] = df['group'].map(group_labels)

plt.rcParams["axes.labelsize"] = 18
g = sns.PairGrid(df, hue="group" , palette=['#2986cd', '#ff3939'])
g.map_diag(sns.histplot, kde = True, element="step", stat="density")
g.map_offdiag(sns.scatterplot, s=150, alpha=0.85)
g.add_legend(title= "")
plt.subplots_adjust(bottom=0.2, left=0.2)
g.fig.supxlabel('NPX, Normalized Protein eXpression', x=0.5, y=0.038, size = 14)
g.fig.supylabel('NPX, Normalized Protein eXpression', x=0.048, y=0.5, size = 14)
plt.savefig('Figure1A_R1.tif', dpi=900)
plt.show()
end


**# 3D Scatter plot
cd "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi"
python

import matplotlib.pyplot as plt 
# change font
plt.rcParams['font.family'] = "Arial"
import numpy as np
import pandas as pd 
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D


df = pd.read_csv("IgANcytochines.csv")
df.dropna(inplace=True)

df.agg(
    {
        "CX3CL1": ["min", "max", "median", "mean"],
		"IL6": ["min", "max", "median", "mean"],
		"IL12B": ["min", "max", "median", "mean"],
	 }	
)


X = df['CX3CL1']
Y = df['IL12B']
Z = df['IL6']

group = df['group']


# Plot 3D plot
fig = plt.figure(figsize=(14,14))
ax =  fig.add_subplot(projection='3d')

# ax.scatter(X, Y, Z, c = group)


ax.scatter(X[group == 0], Y[group == 0], Z[group == 0], c='#2986cd', edgecolor = 'white', linewidth=1.8, marker='o',  s = 1000, alpha = 0.85, label='Controls')
ax.scatter(X[group == 1], Y[group == 1], Z[group == 1], c='#ff3939', edgecolor = 'white', linewidth=1.8,  marker='o', s = 1000, alpha = 0.85, label='IgA Nephropathy')

#legend
ax.legend(prop={'size':25}, loc='upper left')


ax.set_xlabel("C-X3-C motif chemokine ligand 1 (CXC3CL1), NPX", fontsize=16, labelpad = 5)
ax.set_ylabel("Interleukin-12 subunit beta (IL-12B), NPX", fontsize=16)
ax.set_zlabel("Interleukin-6 (IL6), NPX", fontsize=14)

# Optional details
#For x-axis limit
ax.set_xlim(3, 7)
#For y-axis limit
ax.set_ylim(4, 10)
#For z-axis limit
ax.set_zlim(7,1)


ax.tick_params(axis='x', labelsize=14, pad=0.5)
ax.tick_params(axis='y', labelsize=14, pad=1.3)
ax.tick_params(axis='z', labelsize=14, pad=1)


ax.view_init(-140, 60)


# savefig specifies the DPI for the saved figure 
plt.savefig('Figure1B_R1.tif', dpi=900)
plt.show()

end


**# fit multivariable logistic model with selected vars 
cd "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi"
use IgANcytochines, clear
summ $selected_var

* center the variables	
foreach var of varlist $selected_var  {
	 qui summ  `var'
	 qui replace `var' = (`var' - r(mean)) / r(sd)
	 }


di "$selected_var"


**# sensitivity specificity with interactions
logit group c.IL6##c.IL12B##c.CX3CL1
lstat
cap drop pr
predict pr, pr
roctab group pr
local sroc = string(`r(area)', "%4.3f")
local slb = string(`r(lb)', "%4.3f")
local sub = string(`r(ub)', "%4.3f")
lroc, aspectratio(1) note("ROC area = `sroc' (95%CI: `slb' to `sub')")
cap graph export roc_3vars_int.png, replace
cap graph export roc_3vars_int.pdf, replace

**# sensitivity specificity with no interactions
logit group $selected_var
lstat
cap drop pr
predict pr, pr
roctab group pr
local sroc = string(`r(area)', "%4.3f")
local slb = string(`r(lb)', "%4.3f")
local sub = string(`r(ub)', "%4.3f")
lroc, aspectratio(1) note("ROC area = `sroc' (95%CI: `slb' to `sub')")
cap graph export roc_3vars.png, replace
cap graph export roc_3vars.pdf, replace

**# 10-fold cross-validatated AUCROC
cvauroc group $selected_var, seed(123)

**# check 95%CI and P values
qui logit group $selected_var
etable,                             ///
        cstat(_r_b, nformat(%4.2f))                       ///
        cstat(_r_ci, cidelimiter(,) nformat(%6.2f))       ///
        showstars showstarsnote                           /// 
        stars(.05 "*" .01 "**" .001 "***", attach(_r_b))  ///
        mstat(N) mstat(aic) mstat(bic)                    ///
        mstat(pseudo_r2 = e(r2_p))
lincom IL6, or cformat(%3.2f) pformat(%4.3f) sformat(%3.2f)
lincom IL12B, or cformat(%3.2f) pformat(%4.3f) sformat(%3.2f)
lincom CX3CL1, or cformat(%3.2f) pformat(%4.3f) sformat(%3.2f)

**# save estimates and plot ROC curve
foreach var of varlist IL6 IL12B CX3CL1 {
	 qui logit group `var'
	 cap drop pr_`var'
	 predict pr_`var', pr
	 label var pr_`var' "`var'"
	 qui roctab group pr_`var'
	 local sroc_`var' = string(`r(area)', "%4.3f")
	 local slb_`var' = string(`r(lb)', "%4.3f")
	 local sub_`var' = string(`r(ub)', "%4.3f")
	 }
qui logit group IL6 IL12B CX3CL1
	 cap drop pr_ALL3
	 predict pr_ALL3, pr
	 label var pr_ALL3 "IL6+IL12B+CX3CL1"
	 qui roctab group pr_ALL3
	 local sroc_ALL3 = string(`r(area)', "%4.3f")
	 local slb_ALL3 = string(`r(lb)', "%4.3f")
	 local sub_ALL3 = string(`r(ub)', "%4.3f")
	 
roccomp group pr_IL6 pr_IL12B pr_CX3CL1 pr_ALL3, graph ///
plot1opts(lcolor(stc1) mcolor(stc1)) ///
plot2opts(lcolor(stc2) mcolor(stc2)) ///
plot3opts(lcolor(stc3) mcolor(stc3)) ///
plot4opts(lcolor(black) mcolor(black)) ///
legend(order(1 "IL6, ROC area = `sroc_IL6' " 2 "IL12B, ROC area = `sroc_IL12B'" 3 "CX3CL1, ROC area = `sroc_CX3CL1'" 4 "ALL, ROC area = `sroc_ALL3'" ))
cap graph export Figure2.png, replace
cap graph export Figure2.pdf, replace 

**# clinical and lab correlation with cytokines in the multivariable model

**# Matrix graph ggally (PairGrid from seaborn) in Python selected vars /clincal
preserve
clear
import excel "C:\Documenti\Sara Albrandi\Olink_IgAN_clean.xlsx", sheet("Olink + baseline") firstrow
cd "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi"
keep cx3cl1 il12b il6 egfr serum_alb proteinuria
order cx3cl1 il12b il6 egfr serum_alb proteinuria
label var cx3cl1 "CX3CL1, NPX"
label var il12b "IL-12β, NPX"
label var il6 "IL-6, NPX"
label var egfr "eGFR, ml/min/1.73m{sup:2}"
label var serum_alb "Serum Albumin, g/dl"
label var proteinuria "Proteinuria, g/day"
graph matrix cx3cl1 il12b il6 egfr serum_alb proteinuria, half
export excel _all using "clincial_correlations", firstrow(variables) nolabel replace
export delimited _all using clinical_correlations, nolabel quote replace
restore


cd "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi"
python:
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
# change font
plt.rcParams['font.family'] = "Arial"
df = pd.read_csv("clinical_correlations.csv")
df.dropna(inplace=True)
df = df.rename(columns = {"cx3cl1": "CX3CL1, NPX", "il12b": "IL-12B, NPX", "il6": "IL-6, NPX", "egfr": "eGFR, ml/min/1.73m$^{2}$", "serum_alb": "sAlbumin, g/dl", "proteinuria": "Proteinuria, g/day"}) 

plt.rcParams["axes.labelsize"] = 13

g = sns.PairGrid(df, diag_sharey=False, corner=True)
g.map_diag(sns.histplot, color = '#ff3939', alpha = 0.50, bins = 12, kde = True, element="step", stat="density")
g.map_offdiag(sns.scatterplot, s = 150, c='#ff3939', edgecolor = 'white', linewidth=1.8, marker='o', alpha = 0.85)

# g = sns.PairGrid(df, hue="group" , palette=['#2986cd', '#ff3939'])
# g.map_diag(sns.histplot, kde = True, element="step", stat="density")
# g.map_offdiag(sns.scatterplot, s=150, alpha=0.85)
# g.add_legend(title= "")
# plt.subplots_adjust(bottom=0.2, left=0.2)
# g.fig.supxlabel('NPX, Normalized Protein eXpression', x=0.5, y=0.038, size = 14)
# g.fig.supylabel('NPX, Normalized Protein eXpression', x=0.048, y=0.5, size = 14)
plt.savefig('Figure_new.tif', dpi=900)
plt.show()
end


**# compute associated matrix of Spearman's rak correlation and P values  
preserve
clear
import delimited "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi\clinical_correlations.csv"
*tab group
* keep if group == 1
order cx3cl1 il12b il6 egfr serum_alb proteinuria
local vars "cx3cl1 il12b il6 egfr serum_alb proteinuria"
cap drop nmiss
egen nmiss = rowmiss(`vars')
local n_vars: word count `vars'
matrix rho = J(`n_vars',`n_vars',-99)
matrix pval = J(`n_vars',`n_vars',-99)
matrix npts = J(`n_vars',`n_vars',-99)
local i=0
foreach v in `vars' {
           local i= `i' +1
           local j=0
           foreach w in `vars'  {
                   local j= `j'+1
           qui spearman `v' `w' if nmiss==0
           matrix rho[`i',`j'] = r(rho)
		   matrix pval[`i',`j'] = r(p)
		   matrix npts[`i',`j'] = r(N)	   
		   matrix rho[`i',`i'] = .
		   matrix pval[`i',`i'] = .
		   matrix npts[`i',`i'] = .
           }
   }

unab vars : `vars'
matrix rownames rho= `vars'
matrix colnames rho= `vars'
matrix rownames pval= `vars'
matrix colnames pval= `vars'
matrix rownames npts= `vars'
matrix colnames npts= `vars'
matrix list rho,  format(%4.3f)
matrix list pval, format(%4.3f)




clear
import excel "C:\Documenti\Sara Albrandi\Olink_IgAN_clean.xlsx", sheet("Olink + baseline") firstrow

* Spearman's rank correlation
keep if group == 1
foreach clin_var of varlist age gender race proteinuria serum_creat egfr ///
	serum_alb hematuria m s e t c is raas  {
		foreach cytochine of varlist cx3cl1 il12b il6  {
			qui spearman `clin_var' `cytochine'
			di _newline(3) _col(8) "-----> correlation between  `clin_var' and  `cytochine':  rho = " %4.3f r(rho) " P= " %4.3f r(p)
			}
		 }

* gamma regression (univariable)
foreach iga of varlist  iga_elisa gdiga_elisa iggaa_company2  {
	foreach var of varlist cx3cl1 il12b il6  {
			qui summ `var' 
			replace  `var' = (`var' - r(mean)) / r(sd) 
			qui glm `iga' `var' if group == 1, family(gamma) link(log) vce(robust)
			test `var'
			di _newline(3) _col(8) "-----> relation between  `var' and  `iga':  P= " %4.3f r(p)
			}
		 }

		 
* gamma regression (multivariable)
foreach iga of varlist  iga_elisa gdiga_elisa iggaa_company2  {	
						di _newline(3) _col(8) "-----> multivariable relation with  `iga':
						glm `iga' cx3cl1 il12b il6 if group == 1, family(gamma)  link(log) vce(robust)
		 }
		 
		 

* make the plot from multivariable gamma regression
clear
import excel "C:\Documenti\Sara Albrandi\Olink_IgAN_clean.xlsx", sheet("Olink + baseline") firstrow
keep if group == 1
glm iggaa_company2 cx3cl1 il12b il6 if group == 1, family(gamma) link(log) vce(robust)
qui test il12b
local pval = string(r(p), "%4.3f")
margins, at(il12b = (4(0.2)9.8))
marginsplot,  recastci(rarea) ciopts(color(black%10) lcolor(white)) ///
			plotopts(msymbol(i) lcolor(black)) xlabel(4(1)10) ///
			addplot(scatter iggaa_company2 il12b, mcolor(stc2%85) ///
			msize(*1.8) mlcolor(white) mlwidth(*0.7) xlabel(4(1)10)) ///
			legend(off) title(" ")  ///
			ytitle("IgG autoantibodies against gd-IgA1, {&mu}g/ml") ///
			xtitle("IL12B, NPX")  xscale(titlegap(3)) ///
			text(1.8 5 "P = `pval'") name(iga_autoantibodies_il12b, replace)
			
* gd-IgA1 
glm gdiga_elisa cx3cl1 il12b il6 if group == 1, family(gamma) link(log) vce(robust)
qui test il12b
local pval = string(r(p), "%4.3f")
margins, at(il12b = (4(0.2)9.8))
marginsplot,  recastci(rarea) ciopts(color(black%10) lcolor(white)) ///
			plotopts(msymbol(i) lcolor(black)) xlabel(4(1)10) ///
			addplot(scatter gdiga_elisa il12b, mcolor(stc2%85) msize(*1.8) ///
			mlcolor(white) mlwidth(*0.7) xlabel(4(1)10)) ///
			legend(off) title(" ")  ///
			ytitle("gd-IgA1, {&mu}g/ml") ///
			xtitle("IL12B, NPX")  xscale(titlegap(3)) ///
			text(18.5 5 "P = `pval'") name(iga_il12b, replace)
			
graph combine iga_il12b iga_autoantibodies_il12b, ///
name("firstset", replace) xcommon cols(2) 
graph export Figure3.png, replace
graph export Figure3.pdf, replace
			
**# additional correlations 
* Iga, gdIgA, and IgG anti gd-IgA between them and with screat, proteinuria, and eGFR
* foreach var of varlist 

**# compute matrix of Spearman's rak correlation and P values
preserve
tab group
keep if group == 1
order  proteinuria serum_creat egfr iga_elisa gdiga_elisa iggaa_company2
rename proteinuria Proteinuria 
rename serum_creat SCreat 
rename egfr eGFR 
rename iga_elisa IgA 
rename gdiga_elisa gdIgA1 
rename iggaa_company2 IgGantigdIgA1
local vars "Proteinuria SCreat eGFR IgA gdIgA1 IgGantigdIgA1"
cap drop nmiss
egen nmiss = rowmiss(`vars')
local n_vars: word count `vars'
matrix rho = J(`n_vars',`n_vars',-99)
matrix pval = J(`n_vars',`n_vars',-99)
matrix npts = J(`n_vars',`n_vars',-99)
local i=0
foreach v in `vars' {
           local i= `i' +1
           local j=0
           foreach w in `vars'  {
                   local j= `j'+1
           qui spearman `v' `w' if nmiss==0
           matrix rho[`i',`j'] = r(rho)
		   matrix pval[`i',`j'] = r(p)
		   matrix npts[`i',`j'] = r(N)	   
		   matrix rho[`i',`i'] = .
		   matrix pval[`i',`i'] = .
		   matrix npts[`i',`i'] = .
           }
   }

unab vars : `vars'
matrix rownames rho= `vars'
matrix colnames rho= `vars'
matrix rownames pval= `vars'
matrix colnames pval= `vars'
matrix rownames npts= `vars'
matrix colnames npts= `vars'
matrix list rho,  format(%4.3f)
matrix list pval, format(%4.3f)

graph matrix Proteinuria SCreat eGFR IgA gdIgA1 IgGantigdIgA1, half
keep Proteinuria SCreat eGFR IgA gdIgA1 IgGantigdIgA1
export delimited _all using other_correlations, replace

restore
clear
import delimited "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi\other_correlations.csv

**# Matrix graph ggally (PairGrid from seaborn) in Python
cd "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi"
python:
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
# change font
plt.rcParams['font.family'] = "Arial"
df = pd.read_csv("other_correlations.csv")
df.dropna(inplace=True)
df = df.rename(columns = {"Proteinuria": "UProt, g/day", "SCreat": "SCreat, mg/dL", "eGFR": "eGFR, ml/min/1.73m$^{2}$", "IgA": "IgA, ug/ml", "gdIgA1": "gd-IgA1, ug/ml", "IgGantigdIgA1": "IgG anti-gd-IgA1, ug/ml"}) 

plt.rcParams["axes.labelsize"] = 13

g = sns.PairGrid(df, diag_sharey=False, corner=True)
g.map_diag(sns.histplot, color = '#ff3939', alpha = 0.50, bins = 12, kde = True, element="step", stat="density")
g.map_offdiag(sns.scatterplot, s = 150, c='#ff3939', edgecolor = 'white', linewidth=1.8, marker='o', alpha = 0.85)

# g = sns.PairGrid(df, hue="group" , palette=['#2986cd', '#ff3939'])
# g.map_diag(sns.histplot, kde = True, element="step", stat="density")
# g.map_offdiag(sns.scatterplot, s=150, alpha=0.85)
# g.add_legend(title= "")
# plt.subplots_adjust(bottom=0.2, left=0.2)
# g.fig.supxlabel('NPX, Normalized Protein eXpression', x=0.5, y=0.038, size = 14)
# g.fig.supylabel('NPX, Normalized Protein eXpression', x=0.048, y=0.5, size = 14)
plt.savefig('Figure4_R1.tif', dpi=900)
plt.show()
end


**# Matrix graph ggally (PairGrid from seaborn) in Python
cd "C:\Documenti\Alibrandi\IgA ELISA ZORN\Cytochines Cravedi"
python:
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
# change font
plt.rcParams['font.family'] = "Arial"
df = pd.read_csv("other_correlations.csv")
df.dropna(inplace=True)
df = df.rename(columns = {"Proteinuria": "UProt, g/day", "SCreat": "SCreat, mg/dL", "eGFR": "eGFR, ml/min/1.73m$^{2}$", "IgA": "IgA, ug/ml", "gdIgA1": "gd-IgA1, ug/ml", "IgGantigdIgA1": "IgG anti-gd-IgA1, ug/ml"}) 

plt.rcParams["axes.labelsize"] = 13
"""
g = sns.jointplot(df, x='IgG anti-gd-IgA1, ug/ml', y ='eGFR, ml/min/1.73m$^{2}$', kind = 'reg', scatter_kws= {'color':'#ff3939', 'edgecolor':'white', 'alpha' : 0.50, 's': 80}, line_kws= {'color':'#ff3939'}, marginal_kws= {'color':'#ff3939', 'bins': 10})
g.ax_joint.set_yticks([0,15,30,45,60,90,120,150])

plt.savefig('Figure_eGFR_anti-gd-IgA1_R1.tif', dpi=900)
plt.show()
"""

g = sns.jointplot(df, x='IgA, ug/ml', y ='eGFR, ml/min/1.73m$^{2}$', kind = 'reg', scatter_kws= {'color':'#ff3939', 'edgecolor':'white', 'alpha' : 0.50, 's': 80}, line_kws= {'color':'#ff3939'}, marginal_kws= {'color':'#ff3939', 'bins': 10})
g.ax_joint.set_yticks([0,15,30,45,60,90,120,150])
g.ax_joint.grid(True, lw = 0.9, ls = '--', c = '.75', alpha = 0.5)
g.ax_joint.set_xlabel(r'IgA $\mu \$g/ml')
g.ax_joint.text(9000, 95, r'$\rho \$ = -0.339; P=0.020', fontsize = 12) 


plt.savefig('Figure_eGFR_IgA_R1.tif', dpi=300)
plt.show()

"""
g = sns.jointplot(df, x='IgG anti-gd-IgA1, ug/ml', y ='IgA, ug/ml', kind = 'reg', scatter_kws= {'color':'#ff3939', 'edgecolor':'white', 'alpha' : 0.50, 's': 80}, line_kws= {'color':'#ff3939'}, marginal_kws= {'color':'#ff3939', 'bins': 10})


plt.savefig('Figure_anti-gd-IgA1_IgA_R1.tif', dpi=900)
plt.show()
"""

end

**# end analysis

stop

