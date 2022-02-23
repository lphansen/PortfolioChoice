# PortfolioChoice

This repository provides illustration for section 5 and appendix B of papar *Asset Pricing under Smooth Ambiguity in Continuous Time* by [Lars Peter Hansen][id1] [Jianjun Miao][id2].
The lastest version of this paper can be accessed through [this link](#link to be filled).

[id1]: https://larspeterhansen.org/
[id2]: https://people.bu.edu/miaoj/


## Table of contents
- [Prerequisite](#prerequisite)
- [Acessing our project](#acessing)
- [Quick user's guide](#quick-guide)

## <a name="prerequisite"></a>Prerequisite
`Python == 3.8.x`, package manager such as `pip`,  and `jupyter notebook` to checkout the notebooks. 

To go to the documentation of this project, please go to: (link avaliable when published).

This project has been tested in an environment with
> `Python == 3.8.8` and  `jupyter notebook == 6.4.5`

## <a name="acessing"></a>Acessing our project

I. To try the jupyter notebook on Google Colab:
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://github.com/SuriChen/PortfolioChoce.git/PortfolioChoice.ipynb)

II. To store the notebook as well as codes in your local machine. You can do this by following steps below:

1.  Open a Windows command prompt, Mac terminal or Linux terminal and change into the folder you would like to store the files.
 	-  You can do this using the command `cd` in the command prompt (on Windows) or terminal (on Mac and Linux).
        - For example, running `cd 'C:\Users\username\python'` (don’t forget '' around the path name to use an absolute path) would lead me to my designated folder.
     
    ```bash
        cd [folder path name]
    ```

2.  Clone the github repository for the paper
    - If you don’t have github installed, try installing it from this page: https://git-scm.com/download.
    - You can do this by running below in the command prompt:
    
    ```bash
        git clone https://github.com/SuriChen1028/PortfolioChoice.git
    ```
    
3.  Change directories into the ‘Wrestling’ folder and install the required packages for the current user or your initiated virtual environment:
    
    ```bash
        cd PortfolioChoice
        pip install -r requirements.txt
    ```
4. Access the notebooks, run the following under the folder `PortfolioChoice/`:
    
    ```bash
        jupyter notebook
    ```