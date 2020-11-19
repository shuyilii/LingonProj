---
title: "LingonProj Gene Ontology Analysis"
author:
   name: "Shuyi Li"
   email: shuyi.li@med.lu.se
   affiliation: LUDC Bioinformatics Unit
date: "19 November, 2020"
output:
  html_document:
    keep_md: true
---



## Seperate genes based on different patterns

Input for run_select_pattern.sh: RSEM result (GeneMat_HFD_Lingon_LDF.Ebseqresults/GeneMat_HFD_Lingon_LDF.Ebseqresults_FDR_0.05.tab)

Output: Gene list with EnsemblID	and GeneSymbol, classified as diffent patterns
(Output folder: Matrix/pattern)


```bash
#seperate genes based on different patterns
scr/run_select_pattern.sh
#full list of expressed genes
cut -f 1 Matrix/GeneMat_HFD_Lingon_LDF.Ebseqresults | sed $'s/_/\t/' | sed $'s/"//g'|tail -n +2 > Matrix/full_gene_list.csv
```

## GO analysis

### Load data

```r
pattern2.FDR0.05 <- read.table("Matrix/pattern/pattern2_FDR0.05.csv", header = TRUE)
pattern3.FDR0.05 <- read.table("Matrix/pattern/pattern3_FDR0.05.csv", header = TRUE)
pattern4.FDR0.05 <- read.table("Matrix/pattern/pattern4_FDR0.05.csv", header = TRUE)
full_gene <- read.table("Matrix/full_gene_list.csv", header = FALSE)
names(full_gene) <- c("EnsemblID","GeneSymbol")
```


```r
#load packages
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(gprofiler2)
```

### RSEM result summary
Summarize the number of genes classified to each pattern (FDR cutoff = 0.05)

```r
pattern_table <- read.table("Matrix/GeneMat_HFD_Lingon_LDF.Ebseqresults.pattern",header = TRUE)
names(pattern_table) <- c("HDF","Lingon","LFD")
RSEM_summary <-data.frame(c(0,155,40,123,0,318),
                 row.names = c("pattern1.FDR0.05","pattern2.FDR0.05","pattern3.FDR0.05","pattern4.FDR0.05","pattern5.FDR0.05","Sum"))
names(RSEM_summary) <- c("number of genes")
knitr::kable(pattern_table)
```



|         | HDF| Lingon| LFD|
|:--------|---:|------:|---:|
|Pattern1 |   1|      1|   1|
|Pattern2 |   1|      1|   2|
|Pattern3 |   1|      2|   1|
|Pattern4 |   1|      2|   2|
|Pattern5 |   1|      2|   3|

```r
knitr::kable(RSEM_summary)
```



|                 | number of genes|
|:----------------|---------------:|
|pattern1.FDR0.05 |               0|
|pattern2.FDR0.05 |             155|
|pattern3.FDR0.05 |              40|
|pattern4.FDR0.05 |             123|
|pattern5.FDR0.05 |               0|
|Sum              |             318|
For example:

Genes belonging to Pattern 1 equals no differential change between any of the experiments

Genes belonging to Pattern 2 are differentially expressed in LFD compared to HFD & Lingon

Genes in Pattern 3 are DE in Lingon compared to the other experiments

Genes in Pattern 4 have similar expression in Lingon and LFD but are differentially expressed in HFD

### GO analysis using clusterProfiler
(clusterProfiler v3.16.1  For help: https://guangchuangyu.github.io/software/clusterProfiler)

```r
ego<- function(pattern_data){
   ego_pattern <- enrichGO(gene = pattern_data[,1], 
                           universe = full_gene[,1],
                           OrgDb = org.Mm.eg.db, 
                           keyType = "ENSEMBL", 
                           ont = 'ALL',       
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)

   str1 <- "GO_results/clusterProfiler/"
   str2 <- "_GO_enrich.csv"
   file_name <- paste(str1, deparse(substitute(pattern_data)), str2, sep = "")
   write.csv(as.data.frame(ego_pattern),file_name,row.names = F)
   return(ego_pattern)
}
```

#### GO analysis for gene assigned to pattern2.FDR0.05

```r
ego_pattern2.FDR0.05 <- ego(pattern2.FDR0.05)
barplot(ego_pattern2.FDR0.05,showCategory=20,drop=T)
```

![](GO_analysis_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
dotplot(ego_pattern2.FDR0.05,showCategory=20)
```

![](GO_analysis_files/figure-html/unnamed-chunk-6-2.png)<!-- -->


```r
heatplot(ego_pattern2.FDR0.05)
```

![](GO_analysis_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

#### GO analysis for gene assigned to pattern3.FDR0.05

```r
ego_pattern3.FDR0.05 <- ego(pattern3.FDR0.05)
barplot(ego_pattern3.FDR0.05,showCategory=20,drop=T)
```

![](GO_analysis_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
dotplot(ego_pattern3.FDR0.05,showCategory=20)
```

![](GO_analysis_files/figure-html/unnamed-chunk-8-2.png)<!-- -->


```r
heatplot(ego_pattern3.FDR0.05)
```

![](GO_analysis_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

#### GO analysis for gene assigned to pattern4.FDR0.05

```r
ego_pattern4.FDR0.05 <- ego(pattern4.FDR0.05)
barplot(ego_pattern4.FDR0.05,showCategory=20,drop=T)
```

![](GO_analysis_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
dotplot(ego_pattern4.FDR0.05,showCategory=20)
```

![](GO_analysis_files/figure-html/unnamed-chunk-10-2.png)<!-- -->


```r
heatplot(ego_pattern4.FDR0.05)
```

![](GO_analysis_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

### GO analysis using gprofiler
(gprofiler2 version 0.2.0 For help:https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html)

#### GO analysis for gene assigned to pattern2.FDR0.05

```r
gost_pattern2 <- gost(pattern2.FDR0.05[,1], organism ='mmusculus', sources = 'GO', significant = TRUE)
gostplot(gost_pattern2)
```

<!--html_preserve--><div id="htmlwidget-3b3d63488750acf71b85" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-3b3d63488750acf71b85">{"x":{"data":[{"x":[152.704842131299,125.49111906501,137.410448177414,108.559476731527,171.346141863624,139.470082502609,137.42653907058,108.551431284945,171.342119140333,139.466059779317,149.361959076149,179.013452457024,196.347367119646,182.472994487624,161.900787575428,116.427923489497,97.4929649568992,105.582661495895,157.045360562714,140.88608110118,177.255522378684,149.345868182984,128.705274974835,137.406425454123,146.316757544563,111.266769506636,179.009429733732,96.7447384246996,103.458663598038,91.241652962071,186.06528638684,140.918262887511,91.2537211319451,139.462037056026,97.4969876801906,96.4711932408848,159.261881096273,177.21529514577,94.3109908334055,116.419878042914,110.663361012927,91.5393344856342],"y":[5.58223084551869,5.433763947274,4.83219871331951,4.46366380878285,4.39439641132312,4.326784246644,4.03157270062041,3.89493532786194,3.85279855534705,3.77020649902961,3.14572932577151,2.62868499895488,2.55184539970724,2.52409571833704,2.48817794487355,2.48339022832007,2.26579325770963,2.26030241907412,2.15531087641797,2.14532612543389,2.0218729406255,1.9585573160625,1.81562586128092,1.74773281354925,1.66779789973532,1.63772827610143,1.55963499729745,1.52888085237328,1.51742413993923,1.47262331680336,1.46979153502038,1.45729056537986,1.43575501347101,1.41750731373862,1.40663954152605,1.40405686533562,1.36688851837634,1.35126119022025,1.34712776024251,1.33635051928339,1.33514211918594,1.30238272832674],"text":["GO:0065008 (3916) <br> regulation of biological quality <br> 2.617e-06","GO:0042730 (19) <br> fibrinolysis <br> 3.683e-06","GO:0048519 (5541) <br> negative regulation of biological process <br> 1.472e-05","GO:0030195 (45) <br> negative regulation of blood coagulation <br> 3.438e-05","GO:1900047 (46) <br> negative regulation of hemostasis <br> 4.033e-05","GO:0050819 (47) <br> negative regulation of coagulation <br> 4.712e-05","GO:0048523 (5029) <br> negative regulation of cellular process <br> 9.299e-05","GO:0030193 (81) <br> regulation of blood coagulation <br> 1.274e-04","GO:1900046 (82) <br> regulation of hemostasis <br> 1.403e-04","GO:0050818 (84) <br> regulation of coagulation <br> 1.697e-04","GO:0061045 (69) <br> negative regulation of wound healing <br> 7.149e-04","GO:1902042 (29) <br> negative regulation of extrinsic apoptotic signaling pathway via death domain receptors <br> 2.351e-03","GO:2000352 (30) <br> negative regulation of endothelial cell apoptotic process <br> 2.806e-03","GO:1903035 (85) <br> negative regulation of response to wounding <br> 2.992e-03","GO:0080090 (6041) <br> regulation of primary metabolic process <br> 3.250e-03","GO:0034116 (14) <br> positive regulation of heterotypic cell-cell adhesion <br> 3.286e-03","GO:0010883 (60) <br> regulation of lipid storage <br> 5.423e-03","GO:0019915 (93) <br> lipid storage <br> 5.492e-03","GO:0071704 (11084) <br> organic substance metabolic process <br> 6.993e-03","GO:0051235 (347) <br> maintenance of location <br> 7.156e-03","GO:1901575 (1967) <br> organic substance catabolic process <br> 9.509e-03","GO:0061041 (145) <br> regulation of wound healing <br> 1.100e-02","GO:0044262 (312) <br> cellular carbohydrate metabolic process <br> 1.529e-02","GO:0048518 (6323) <br> positive regulation of biological process <br> 1.788e-02","GO:0060255 (6489) <br> regulation of macromolecule metabolic process <br> 2.149e-02","GO:0031639 (22) <br> plasminogen activation <br> 2.303e-02","GO:1902041 (47) <br> regulation of extrinsic apoptotic signaling pathway via death domain receptors <br> 2.757e-02","GO:0010675 (166) <br> regulation of cellular carbohydrate metabolic process <br> 2.959e-02","GO:0019222 (6958) <br> regulation of metabolic process <br> 3.038e-02","GO:0007596 (169) <br> blood coagulation <br> 3.368e-02","GO:1904036 (49) <br> negative regulation of epithelial cell apoptotic process <br> 3.390e-02","GO:0051246 (2856) <br> regulation of protein metabolic process <br> 3.489e-02","GO:0007599 (171) <br> hemostasis <br> 3.666e-02","GO:0050817 (172) <br> coagulation <br> 3.824e-02","GO:0010884 (25) <br> positive regulation of lipid storage <br> 3.921e-02","GO:0010605 (2984) <br> negative regulation of macromolecule metabolic process <br> 3.944e-02","GO:0072378 (9) <br> blood coagulation, fibrin clot formation <br> 4.296e-02","GO:1901564 (6486) <br> organonitrogen compound metabolic process <br> 4.454e-02","GO:0009894 (880) <br> regulation of catabolic process <br> 4.496e-02","GO:0034114 (26) <br> regulation of heterotypic cell-cell adhesion <br> 4.609e-02","GO:0031323 (6225) <br> regulation of cellular metabolic process <br> 4.622e-02","GO:0008152 (11640) <br> metabolic process <br> 4.984e-02"],"key":["GO:0065008","GO:0042730","GO:0048519","GO:0030195","GO:1900047","GO:0050819","GO:0048523","GO:0030193","GO:1900046","GO:0050818","GO:0061045","GO:1902042","GO:2000352","GO:1903035","GO:0080090","GO:0034116","GO:0010883","GO:0019915","GO:0071704","GO:0051235","GO:1901575","GO:0061041","GO:0044262","GO:0048518","GO:0060255","GO:0031639","GO:1902041","GO:0010675","GO:0019222","GO:0007596","GO:1904036","GO:0051246","GO:0007599","GO:0050817","GO:0010884","GO:0010605","GO:0072378","GO:1901564","GO:0009894","GO:0034114","GO:0031323","GO:0008152"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(255,153,0,1)","opacity":0.8,"size":[20.9825214744859,11.0344885623109,21.4371807931721,13.371498838758,13.4236723767158,13.4744520389306,21.3113712703715,14.6811714456067,14.7068385270612,14.757072105558,14.3400275960802,12.2614639455316,12.3522331566276,14.7816581511172,21.5485378624855,9.99964881006429,14.0334172271397,14.9667594307067,22.3119836815599,17.3912647470435,20.0430276903694,15.8391513561262,17.2117447350625,21.6070639047992,21.6402214245883,11.4820446595472,13.4744520389306,16.0925528213109,21.7292042930973,16.1257193760994,13.5720986475133,20.5583722310935,16.1474565214948,16.1582157927639,11.8520727978321,20.6179242097883,8.08239076487458,21.6396302927545,18.8716767477501,11.9622477857052,21.5870480322216,22.3721842377421],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(255,153,0,1)"}},"hoveron":"points","set":"SharedData1bfe4711","name":"GO:BP","legendgroup":"GO:BP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[56.2489000321589,68.0296870799672,56.39374577455,66.1024340075969,56.8323064945675,56.2448765393147],"y":[2.64881786771709,2.40798941572576,1.97905766569017,1.83351745128742,1.67209133210063,1.37071104862563],"text":["GO:0005577 (6) <br> fibrinogen complex <br> 2.245e-03","GO:0072562 (7) <br> blood microparticle <br> 3.909e-03","GO:0005615 (1875) <br> extracellular space <br> 1.049e-02","GO:0062023 (359) <br> collagen-containing extracellular matrix <br> 1.467e-02","GO:0005737 (11077) <br> cytoplasm <br> 2.128e-02","GO:0005576 (2780) <br> extracellular region <br> 4.259e-02"],"key":["GO:0005577","GO:0072562","GO:0005615","GO:0062023","GO:0005737","GO:0005576"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(16,150,24,1)","opacity":0.8,"size":[3.77952755905512,6.43262693801819,19.9756429262658,17.4481709405018,22.3112053826337,20.5216318026302],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(16,150,24,1)"}},"hoveron":"points","set":"SharedData1bfe4711","name":"GO:CC","legendgroup":"GO:CC","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[9.62745719752277,9.55906754965469,9.57113631104317,10.4360642105513],"y":[4.14651166761475,3.99973806126618,3.09032705180993,1.83129921444156],"text":["GO:0005515 (10383) <br> protein binding <br> 7.137e-05","GO:0005488 (14952) <br> binding <br> 1.001e-04","GO:0005496 (120) <br> steroid binding <br> 8.122e-04","GO:0008289 (798) <br> lipid binding <br> 1.475e-02"],"key":["GO:0005515","GO:0005488","GO:0005496","GO:0008289"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(220,57,18,1)","opacity":0.8,"size":[22.2313204003894,22.6771653543307,15.4753911606817,18.7229725186289],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(220,57,18,1)"}},"hoveron":"points","set":"SharedData1bfe4711","name":"GO:MF","legendgroup":"GO:MF","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[2,51.1962943399293],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(220,57,18,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[55.2188858640445,73.1757344276949],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(16,150,24,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[77.1983259518102,200.004022591524],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(255,153,0,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[180],"y":[16.2],"text":"values above this threshold are capped","hovertext":"","textfont":{"size":7.55905511811024,"color":"rgba(190,190,190,1)"},"type":"scatter","mode":"text","hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,210],"y":[16,16],"text":"","type":"scatter","mode":"lines","line":{"width":0.755905511811024,"color":"rgba(190,190,190,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":29.2835201328352,"r":6.6417600664176,"b":55.8781661388202,"l":61.0377750103778},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0,210],"tickmode":"array","ticktext":["GO:MF","GO:CC","GO:BP"],"tickvals":[26.5981471699647,64.1973101458697,138.601174271667],"categoryorder":"array","categoryarray":["GO:MF","GO:CC","GO:BP"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.65296803652968,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":10.6268161062682},"tickangle":-45,"showline":true,"linecolor":"rgba(190,190,190,1)","linewidth":0.132835201328352,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-1,18.5],"tickmode":"array","ticktext":["0","2","4","6","8","10","12","14",">16"],"tickvals":[0,2,4,6,8,10,12,14,16],"categoryorder":"array","categoryarray":["0","2","4","6","8","10","12","14",">16"],"nticks":null,"ticks":"outside","tickcolor":"rgba(190,190,190,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":true,"linecolor":"rgba(190,190,190,1)","linewidth":0.66417600664176,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"annotations":[{"text":"-log10(p-adj)","x":-0.0414121646998359,"y":0.5,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":13.2835201328352},"xref":"paper","yref":"paper","textangle":-90,"xanchor":"right","yanchor":"center","annotationType":"axis"},{"text":"query_1","x":0.5,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(169,169,169,1)","family":"","size":13.2835201328352},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"}],"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":24.9730178497302,"yanchor":1,"ysizemode":"pixel"}],"showlegend":false,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895}},"hovermode":"closest","barmode":"relative","dragmode":"zoom"},"config":{"doubleClick":"reset","showSendToCloud":false},"source":"A","attrs":{"707f7906a84f":{"colour":{},"size":{},"alpha":{},"key":{},"x":{},"y":{},"text":{},"type":"scatter"},"707f2414189c":{"x":{},"xend":{},"y":{},"yend":{}},"707f235bbc45":{"x":{},"xend":{},"y":{},"yend":{}},"707f5ba56825":{"x":{},"xend":{},"y":{},"yend":{}},"707f48588424":{"x":{},"y":{}},"707f2b4b7c97":{"yintercept":{}}},"cur_data":"707f7906a84f","visdat":{"707f7906a84f":["function (y) ","x"],"707f2414189c":["function (y) ","x"],"707f235bbc45":["function (y) ","x"],"707f5ba56825":["function (y) ","x"],"707f48588424":["function (y) ","x"],"707f2b4b7c97":["function (y) ","x"]},"highlight":{"on":"plotly_click","off":"plotly_doubleclick","persistent":false,"dynamic":false,"color":null,"selectize":false,"defaultValues":null,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0,"ctGroups":["SharedData1bfe4711"]},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

#### GO analysis for gene assigned to pattern3.FDR0.05

```r
gost_pattern3 <- gost(pattern3.FDR0.05[,1], organism ='mmusculus', sources = 'GO', significant = TRUE)
gostplot(gost_pattern3)
```

<!--html_preserve--><div id="htmlwidget-9d92e87a090f32757aac" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-9d92e87a090f32757aac">{"x":{"data":[{"x":[124.835415168513,87.609133829941,129.352933424749,127.269162759807,97.0826471811768],"y":[2.14012551483488,2.02229090533078,2.02229090533078,1.79378190106033,1.37401499878293],"text":["GO:0042438 (22) <br> melanin biosynthetic process <br> 7.242e-03","GO:0006582 (24) <br> melanin metabolic process <br> 9.500e-03","GO:0044550 (24) <br> secondary metabolite biosynthetic process <br> 9.500e-03","GO:0043473 (92) <br> pigmentation <br> 1.608e-02","GO:0010765 (39) <br> positive regulation of sodium ion transport <br> 4.227e-02"],"key":["GO:0042438","GO:0006582","GO:0044550","GO:0043473","GO:0010765"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(255,153,0,1)","opacity":0.8,"size":[18.7350049677963,19.0039285617227,19.0039285617227,22.6771653543307,20.4249070525632],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(255,153,0,1)"}},"hoveron":"points","set":"SharedData3b3d6348","name":"GO:BP","legendgroup":"GO:BP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[6.48153339559091,16.9612412012591],"y":[3.02419278425077,2.5475241339587],"text":["GO:0004503 (2) <br> monophenol monooxygenase activity <br> 9.458e-04","GO:0016716 (3) <br> oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, another compound as one donor, and incorporation of one atom of oxygen <br> 2.834e-03"],"key":["GO:0004503","GO:0016716"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(220,57,18,1)","opacity":0.8,"size":[3.77952755905512,9.92934068462193],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(220,57,18,1)"}},"hoveron":"points","set":"SharedData3b3d6348","name":"GO:MF","legendgroup":"GO:MF","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[2,51.1962943399293],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(220,57,18,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[55.2188858640445,73.1757344276949],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(16,150,24,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[77.1983259518102,200.004022591524],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(255,153,0,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[180],"y":[16.2],"text":"values above this threshold are capped","hovertext":"","textfont":{"size":7.55905511811024,"color":"rgba(190,190,190,1)"},"type":"scatter","mode":"text","hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,210],"y":[16,16],"text":"","type":"scatter","mode":"lines","line":{"width":0.755905511811024,"color":"rgba(190,190,190,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":29.2835201328352,"r":6.6417600664176,"b":55.8781661388202,"l":61.0377750103778},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0,210],"tickmode":"array","ticktext":["GO:MF","GO:CC","GO:BP"],"tickvals":[26.5981471699647,64.1973101458697,138.601174271667],"categoryorder":"array","categoryarray":["GO:MF","GO:CC","GO:BP"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.65296803652968,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":10.6268161062682},"tickangle":-45,"showline":true,"linecolor":"rgba(190,190,190,1)","linewidth":0.132835201328352,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-1,18.5],"tickmode":"array","ticktext":["0","2","4","6","8","10","12","14",">16"],"tickvals":[0,2,4,6,8,10,12,14,16],"categoryorder":"array","categoryarray":["0","2","4","6","8","10","12","14",">16"],"nticks":null,"ticks":"outside","tickcolor":"rgba(190,190,190,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":true,"linecolor":"rgba(190,190,190,1)","linewidth":0.66417600664176,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"annotations":[{"text":"-log10(p-adj)","x":-0.0414121646998359,"y":0.5,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":13.2835201328352},"xref":"paper","yref":"paper","textangle":-90,"xanchor":"right","yanchor":"center","annotationType":"axis"},{"text":"query_1","x":0.5,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(169,169,169,1)","family":"","size":13.2835201328352},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"}],"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":24.9730178497302,"yanchor":1,"ysizemode":"pixel"}],"showlegend":false,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895}},"hovermode":"closest","barmode":"relative","dragmode":"zoom"},"config":{"doubleClick":"reset","showSendToCloud":false},"source":"A","attrs":{"707f4b5f4d85":{"colour":{},"size":{},"alpha":{},"key":{},"x":{},"y":{},"text":{},"type":"scatter"},"707f5dda7d6b":{"x":{},"xend":{},"y":{},"yend":{}},"707f33602bf0":{"x":{},"xend":{},"y":{},"yend":{}},"707f6ee4b3e9":{"x":{},"xend":{},"y":{},"yend":{}},"707f68dfbedf":{"x":{},"y":{}},"707f396c5843":{"yintercept":{}}},"cur_data":"707f4b5f4d85","visdat":{"707f4b5f4d85":["function (y) ","x"],"707f5dda7d6b":["function (y) ","x"],"707f33602bf0":["function (y) ","x"],"707f6ee4b3e9":["function (y) ","x"],"707f68dfbedf":["function (y) ","x"],"707f396c5843":["function (y) ","x"]},"highlight":{"on":"plotly_click","off":"plotly_doubleclick","persistent":false,"dynamic":false,"color":null,"selectize":false,"defaultValues":null,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0,"ctGroups":["SharedData3b3d6348"]},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

#### GO analysis for gene assigned to pattern4.FDR0.05

```r
gost_pattern4 <- gost(pattern4.FDR0.05[,1], organism ='mmusculus', sources = 'GO', significant = TRUE)
gostplot(gost_pattern4)
```

<!--html_preserve--><div id="htmlwidget-69120273c9f431db580b" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-69120273c9f431db580b">{"x":{"data":[{"x":[142.24173885038,138.697719630661,123.608484564637,88.860200773565,84.9903409672427,146.469621029636,146.437439243305,79.0165968795208,118.209989907585,89.1498368505454,79.1493467481369,108.808885575594,166.277510516466,77.624734620698,100.795620779135,89.1538595738368,129.875887452631,140.990671906756,146.381121117226,107.948022791235,77.206371398393,118.873739250665,90.1233358870631,112.883904269777,129.706933074392,159.201540246902,142.338284209373,138.645424227872,140.705058553067,166.543010253698,89.6245181989301,108.000318194024,89.1297232340885,166.808509990931,137.692038807812,195.52270884491,79.2177330440906,112.887926993069,158.767086131431,158.083223171894],"y":[10.3940392699991,10.3940392699991,8.98626021959794,8.69574275952298,7.70044237860438,6.81781316060113,6.59826525099157,6.59826525099157,6.32735507359577,6.25419587873689,4.93467104505643,4.72502533497339,4.72502533497339,4.22772518618456,3.89647991726821,3.63812071690459,3.32195706584995,3.28911652679177,2.8172127965466,2.56605364782912,2.56138685645222,2.33566423716536,2.2979509155322,2.20635304406003,2.20061949455408,2.16405343360781,2.08886248250516,1.99398399000146,1.81965963975559,1.81121712446427,1.66151085374847,1.60422115982976,1.52894188348894,1.50730859376086,1.49665141259908,1.44978759255888,1.4225948963361,1.41119094285255,1.31099749195991,1.31099749195991],"text":["GO:0051674 (1633) <br> localization of cell <br> 4.036e-11","GO:0048870 (1633) <br> cell motility <br> 4.036e-11","GO:0040011 (1815) <br> locomotion <br> 1.032e-09","GO:0006928 (2060) <br> movement of cell or subcellular component <br> 2.015e-09","GO:0003341 (174) <br> cilium movement <br> 1.993e-08","GO:0060294 (126) <br> cilium movement involved in cell motility <br> 1.521e-07","GO:0060285 (132) <br> cilium-dependent cell motility <br> 2.522e-07","GO:0001539 (132) <br> cilium or flagellum-dependent cell motility <br> 2.522e-07","GO:0035082 (75) <br> axoneme assembly <br> 4.706e-07","GO:0007017 (882) <br> microtubule-based process <br> 5.569e-07","GO:0001578 (107) <br> microtubule bundle formation <br> 1.162e-05","GO:0030317 (113) <br> flagellated sperm motility <br> 1.884e-05","GO:0097722 (113) <br> sperm motility <br> 1.884e-05","GO:0000226 (612) <br> microtubule cytoskeleton organization <br> 5.919e-05","GO:0016477 (1479) <br> cell migration <br> 1.269e-04","GO:0007018 (376) <br> microtubule-based movement <br> 2.301e-04","GO:0044782 (334) <br> cilium organization <br> 4.765e-04","GO:0051270 (1075) <br> regulation of cellular component movement <br> 5.139e-04","GO:0060271 (306) <br> cilium assembly <br> 1.523e-03","GO:0022414 (1734) <br> reproductive process <br> 2.716e-03","GO:0000003 (1735) <br> reproduction <br> 2.745e-03","GO:0035295 (1129) <br> tube development <br> 4.617e-03","GO:0007288 (19) <br> sperm axoneme assembly <br> 5.036e-03","GO:0032502 (6516) <br> developmental process <br> 6.218e-03","GO:0044703 (1257) <br> multi-organism reproductive process <br> 6.301e-03","GO:0072359 (1159) <br> circulatory system development <br> 6.854e-03","GO:0051704 (1278) <br> multi-organism process <br> 8.150e-03","GO:0048856 (5958) <br> anatomical structure development <br> 1.014e-02","GO:0051179 (6194) <br> localization <br> 1.515e-02","GO:0098609 (819) <br> cell-cell adhesion <br> 1.544e-02","GO:0007155 (1363) <br> cell adhesion <br> 2.180e-02","GO:0022610 (1375) <br> biological adhesion <br> 2.488e-02","GO:0007010 (1391) <br> cytoskeleton organization <br> 2.958e-02","GO:0098742 (273) <br> cell-cell adhesion via plasma-membrane adhesion molecules <br> 3.110e-02","GO:0048609 (1072) <br> multicellular organismal reproductive process <br> 3.187e-02","GO:2000145 (977) <br> regulation of cell motility <br> 3.550e-02","GO:0001667 (429) <br> ameboidal-type cell migration <br> 3.779e-02","GO:0032504 (1088) <br> multicellular organism reproduction <br> 3.880e-02","GO:0072240 (2) <br> metanephric DCT cell differentiation <br> 4.887e-02","GO:0072069 (2) <br> DCT cell differentiation <br> 4.887e-02"],"key":["GO:0051674","GO:0048870","GO:0040011","GO:0006928","GO:0003341","GO:0060294","GO:0060285","GO:0001539","GO:0035082","GO:0007017","GO:0001578","GO:0030317","GO:0097722","GO:0000226","GO:0016477","GO:0007018","GO:0044782","GO:0051270","GO:0060271","GO:0022414","GO:0000003","GO:0035295","GO:0007288","GO:0032502","GO:0044703","GO:0072359","GO:0051704","GO:0048856","GO:0051179","GO:0098609","GO:0007155","GO:0022610","GO:0007010","GO:0098742","GO:0048609","GO:2000145","GO:0001667","GO:0032504","GO:0072240","GO:0072069"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(255,153,0,1)","opacity":0.8,"size":[20.5098153897692,20.5098153897692,20.6411290248971,20.7971478311373,17.4334691870896,16.9307970719294,17.0044236288471,17.0044236288471,16.0798833020489,19.7228099104495,16.6687677408511,16.7568187588049,16.7568187588049,19.2369552054154,20.3857782310248,18.564526483353,18.3963519153922,19.9798037849132,18.2707841695982,20.5845192314384,20.5852352997547,20.0428282284777,13.4738853588777,22.1553609846659,20.1800823905409,20.0764514273406,20.2011550868095,22.0533881871558,22.097705002862,19.6254924009223,20.2827968585479,20.2938794442452,20.3084953345489,18.1054773718101,19.9762027971679,19.856173916041,18.7495321260677,19.9952836973811,3.77952755905512,3.77952755905512],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(255,153,0,1)"}},"hoveron":"points","set":"SharedData8750acf7","name":"GO:BP","legendgroup":"GO:BP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[70.8783200136585,63.7205262438324,57.5122767852367,68.3716839717239,57.5082532923925,57.2628202288965,70.5926520217206,60.7914234532571,61.5760045578755,69.4338860825919,58.6589722458328,58.7595595669377,65.1287487393013,63.8814659576003,63.8653719862235,59.6769159354146,57.3151256358711,69.763812495816,70.089715416196,59.419412393386,72.8739724643802],"y":[9.84954394409042,9.29835065635966,7.53872251479426,7.46433726811301,6.73183845855917,5.70606426133768,5.04452387721996,5.03992717779699,4.89865843079297,3.54549506070333,2.81444025436855,2.59727729267921,2.49910703862161,2.46282104288875,2.42648376963498,2.25621072033663,2.06862045310679,1.63097296653717,1.56595829895207,1.5653971940524,1.30153509689639],"text":["GO:0120025 (2397) <br> plasma membrane bounded cell projection <br> 1.414e-10","GO:0042995 (2488) <br> cell projection <br> 5.031e-10","GO:0005930 (126) <br> axoneme <br> 2.893e-08","GO:0097014 (128) <br> ciliary plasm <br> 3.433e-08","GO:0005929 (669) <br> cilium <br> 1.854e-07","GO:0005856 (2247) <br> cytoskeleton <br> 1.968e-06","GO:0099568 (270) <br> cytoplasmic region <br> 9.026e-06","GO:0031514 (217) <br> motile cilium <br> 9.122e-06","GO:0032838 (224) <br> plasma membrane bounded cell projection cytoplasm <br> 1.263e-05","GO:0097729 (139) <br> 9+2 motile cilium <br> 2.848e-04","GO:0015630 (1272) <br> microtubule cytoskeleton <br> 1.533e-03","GO:0016324 (384) <br> apical plasma membrane <br> 2.528e-03","GO:0045177 (471) <br> apical part of cell <br> 3.169e-03","GO:0043232 (4476) <br> intracellular non-membrane-bounded organelle <br> 3.445e-03","GO:0043228 (4490) <br> non-membrane-bounded organelle <br> 3.746e-03","GO:0030286 (60) <br> dynein complex <br> 5.544e-03","GO:0005874 (438) <br> microtubule <br> 8.538e-03","GO:0098590 (1298) <br> plasma membrane region <br> 2.339e-02","GO:0098862 (193) <br> cluster of actin-based cell projections <br> 2.717e-02","GO:0030054 (2034) <br> cell junction <br> 2.720e-02","GO:1990716 (4) <br> axonemal central apparatus <br> 4.994e-02"],"key":["GO:0120025","GO:0042995","GO:0005930","GO:0097014","GO:0005929","GO:0005856","GO:0099568","GO:0031514","GO:0032838","GO:0097729","GO:0015630","GO:0016324","GO:0045177","GO:0043232","GO:0043228","GO:0030286","GO:0005874","GO:0098590","GO:0098862","GO:0030054","GO:1990716"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(16,150,24,1)","opacity":0.8,"size":[20.9819777840206,21.0271294189765,16.9307970719294,16.9557678323324,19.356740133949,20.9033892301901,18.0893689515135,17.7669938230138,17.8142903819691,17.0857266124811,20.1951726090997,18.5942186281003,18.8791992528986,21.7237369433616,21.727369125967,15.6952140480049,18.7784525115574,20.2208803230326,17.5909857995176,20.7815617208014,9.15870512536166],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(16,150,24,1)"}},"hoveron":"points","set":"SharedData8750acf7","name":"GO:CC","legendgroup":"GO:CC","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[11.2044420189515,39.4815499521728,31.2466317647629,51.0796296465073,9.62745719752277],"y":[2.76355153192932,1.99450340476297,1.76554245321608,1.66212326703883,1.32427218675254],"text":["GO:0008569 (19) <br> ATP-dependent microtubule motor activity, minus-end-directed <br> 1.724e-03","GO:0051959 (29) <br> dynein light intermediate chain binding <br> 1.013e-02","GO:0045505 (33) <br> dynein intermediate chain binding <br> 1.716e-02","GO:1990939 (35) <br> ATP-dependent microtubule motor activity <br> 2.177e-02","GO:0005515 (10383) <br> protein binding <br> 4.739e-02"],"key":["GO:0008569","GO:0051959","GO:0045505","GO:1990939","GO:0005515"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(220,57,18,1)","opacity":0.8,"size":[13.4738853588777,14.3451703203736,14.5974189190701,14.7103592075227,22.6771653543307],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(220,57,18,1)"}},"hoveron":"points","set":"SharedData8750acf7","name":"GO:MF","legendgroup":"GO:MF","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[2,51.1962943399293],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(220,57,18,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[55.2188858640445,73.1757344276949],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(16,150,24,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[77.1983259518102,200.004022591524],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(255,153,0,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[180],"y":[16.2],"text":"values above this threshold are capped","hovertext":"","textfont":{"size":7.55905511811024,"color":"rgba(190,190,190,1)"},"type":"scatter","mode":"text","hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,210],"y":[16,16],"text":"","type":"scatter","mode":"lines","line":{"width":0.755905511811024,"color":"rgba(190,190,190,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":29.2835201328352,"r":6.6417600664176,"b":55.8781661388202,"l":61.0377750103778},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0,210],"tickmode":"array","ticktext":["GO:MF","GO:CC","GO:BP"],"tickvals":[26.5981471699647,64.1973101458697,138.601174271667],"categoryorder":"array","categoryarray":["GO:MF","GO:CC","GO:BP"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.65296803652968,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":10.6268161062682},"tickangle":-45,"showline":true,"linecolor":"rgba(190,190,190,1)","linewidth":0.132835201328352,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-1,18.5],"tickmode":"array","ticktext":["0","2","4","6","8","10","12","14",">16"],"tickvals":[0,2,4,6,8,10,12,14,16],"categoryorder":"array","categoryarray":["0","2","4","6","8","10","12","14",">16"],"nticks":null,"ticks":"outside","tickcolor":"rgba(190,190,190,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":true,"linecolor":"rgba(190,190,190,1)","linewidth":0.66417600664176,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"annotations":[{"text":"-log10(p-adj)","x":-0.0414121646998359,"y":0.5,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":13.2835201328352},"xref":"paper","yref":"paper","textangle":-90,"xanchor":"right","yanchor":"center","annotationType":"axis"},{"text":"query_1","x":0.5,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(169,169,169,1)","family":"","size":13.2835201328352},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"}],"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":24.9730178497302,"yanchor":1,"ysizemode":"pixel"}],"showlegend":false,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895}},"hovermode":"closest","barmode":"relative","dragmode":"zoom"},"config":{"doubleClick":"reset","showSendToCloud":false},"source":"A","attrs":{"707f5c2a4ab9":{"colour":{},"size":{},"alpha":{},"key":{},"x":{},"y":{},"text":{},"type":"scatter"},"707f5c8fe6f4":{"x":{},"xend":{},"y":{},"yend":{}},"707f6b83cca5":{"x":{},"xend":{},"y":{},"yend":{}},"707f19f09bc8":{"x":{},"xend":{},"y":{},"yend":{}},"707f38374c6":{"x":{},"y":{}},"707f2b6774f7":{"yintercept":{}}},"cur_data":"707f5c2a4ab9","visdat":{"707f5c2a4ab9":["function (y) ","x"],"707f5c8fe6f4":["function (y) ","x"],"707f6b83cca5":["function (y) ","x"],"707f19f09bc8":["function (y) ","x"],"707f38374c6":["function (y) ","x"],"707f2b6774f7":["function (y) ","x"]},"highlight":{"on":"plotly_click","off":"plotly_doubleclick","persistent":false,"dynamic":false,"color":null,"selectize":false,"defaultValues":null,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0,"ctGroups":["SharedData8750acf7"]},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

### GO analysis output

#### clusterProfiler

[clusterProfiler GO analysis output](../GO_results/clusterProfiler)

clusterProfiler GO analysis output files are in the ./GO_results/clusterProfiler

#### gprofiler

[gprofiler GO analysis output - pattern2.FDR0.05](https://biit.cs.ut.ee/gplink/l/faZQDvB8SE)

[gprofiler GO analysis output - pattern3.FDR0.05](https://biit.cs.ut.ee/gplink/l/0_2ljZG4SM)

[gprofiler GO analysis output - pattern4.FDR0.05](https://biit.cs.ut.ee/gplink/l/0y1ajpx2Sn)

gprofiler GO analysis output files are in the ./GO_results/gprofiler
