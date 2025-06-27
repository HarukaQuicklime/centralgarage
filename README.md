-2025.6.27 update
-New_Vyazovkin_Ea_Hikari_ver. 1. 1-For isoconversional kinetic analysis

The Excel file is a convenient tool for kinetic analysis.

About this file

This Excel-based tool is designed as a plug-and-play solution for thermal analysis. No coding experience is required — simply paste your data and let the calculations run automatically.

It supports activation energy calculation using the Vyazovkin method with three, four, or five different heating rates. The temperature integral is approximated using the Cai approximation (see references for details).

The activation energy (Ea) is evaluated in the range of 50 to 500 kJ/mol with a resolution of 0.01 kJ/mol.

Results obtained from this file show only minor deviations from other widely used isoconversional methods (FWO, KAS, Friedman, and Starink), indicating high stability and good consistency.

You may use 3 to 5 heating rates for the calculation. However, using 4 or 5 is recommended to reduce numerical error and improve accuracy.

Instructions:
1. Enter the heating rates (in K/min), e.g., 5, 10, 15, 20, 25. If only three rates are available, set the fourth and fifth values to 0.
2. Input the temperature (in K) corresponding to each conversion level (α).
3. The activation energy will be calculated automatically.
4. For any unused conversion level, enter 0 in the corresponding cells.

All editable fields are highlighted in light blue (or in a color distinct from black) to ensure clarity and accessibility. The color scheme is designed to be safe for most types of color vision deficiency.


Developed by **Hikari Quicklime** (nickname of the first author).

If you find this tool helpful in your research, please consider citing the corresponding article. 🥰

**Note:
**This file is relatively large (~70 MB) . In the future, the author (Hikari) plans to develop a more compact, Python-based application with an interactive graphical user interface (GUI), aiming to improve speed, flexibility, and ease of use.

(Yes... another flag has been raised 😂)


you can download the file from https://drive.google.com/file/d/1ss4nfaRVdHHvHJU2lQEUCMoEXKLjU9nS/view?usp=sharing![image](https://github.com/user-attachments/assets/3a76d18c-52b5-4910-a6ad-0a77c0bb648a)


(old)-2022
代码使用Matlab语言编写，MATLAB R2019b运行OK。
这是一个三高斯分布式活化能模型。
你可以用这个模型计算材料热解的动力学参数。
本研究使用的指前因子为10^17.5，升温速率为1/6 oC/s。
初学者编写，轻喷。
