conda create -n plotly python=3

conda activate plotly
pip install plotly
pip install -U kaleido
conda install -c conda-forge pyinstaller
# pip install --upgrade pyinstaller
pyinstaller -F plotly_go.py -p /home/dozeduck/workspace/anaconda3/envs/plotly/lib/python3.11/site-packages

# for R packages
# conda install pandas
conda install rpy2
python
import rpy2.robjects as robjects
# 启动 R 环境
r = robjects.r
# 定义 R 代码字符串
r_code = '''
install.packages("ggplot2")
install.packages("hexbin")
install.packages("KernSmooth")
install.packages("readxl")
'''
# 执行 R 代码
r(r_code)

如果上面的代码不行就尝试下面的：

import rpy2.robjects as ro

# 设置CRAN镜像
ro.r('options(repos = c(CRAN = "http://cran.rstudio.com/"))')

# 安装ggplot2包
ro.r('install.packages("ggplot2")')

