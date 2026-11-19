# 1000 设计的追加筛选与完整 IgG 序列装配复核

完成时间：2026-11-19T03:12:55.498673+00:00。本轮在当前节点完成 CPU 分析，无新增设计生成，无外部提交。

从 1000 个设计（998 个独立 VH/VL 对）出发，原有几何/序列审查池 237 个；追加复核后保留 211 个，选出 20 个优先候选和 20 个备选。这里的“优先”表示值得实验测试，并非已经验证结合。原先包含 OF3 的严格全门槛通过数仍为 0。

## 是否降低标准

- 原有几何、目标糖基化位点净距、CDR 糖基化、额外 Cys 和新颖性阈值未放宽。OF3 继续仅作参考：已知阳性对照未校准通过，不能用作可靠淘汰依据；未测与测过不通过保持不同标签。
- 本轮新增的保守筛选排除 24 个含连续至少 5 个 AVILMFWY 残基的 CDR 候选，以及 2 个 Fv 重轻链间存在 <2 Å 重原子冲突的候选。这些是计算优选条件，不是比赛官方规则。
- 211 个中，36 个的六条 CDR 均不含 N[GS]/D[GS]；优先挑选这一组。施加六 CDR 相似度 <0.8 的多样性要求后，得到 29 个此类候选，足够组成优先 20 个，另有 9 个进入备选。
- 剩余 11 个备选允许 N[GS]/D[GS]，逐位保留风险标记，没有修改序列。这是唯一用于补齐备选的软偏好放宽。即使无这些基序，也不能称为无降解/表达风险；M/W、DP 等其他提示仍列在每条记录中。
- 为优先保留较少基序风险和多样性，不强制四个生成条件各 5 个；优先名单四条件数量为 7、6、2、5。相较旧优先 20，替换 17 个；逐条去向见 `reports/previous_priority_disposition.csv`。

## 新颖性与序列完整性

- 对设计与本地冻结参考库共 4020 条去重链序列重新进行 ANARCII / IMGT 编号，全部编号成功。核对 VH / kappa VL 链型、完整可变区覆盖、框架不变、六 CDR 与设计注释一致、IMGT 23/104 保守 Cys，以及整个可变区的 N-X-S/T 基序（X != P）。
- 扩展比较本地 TheraSAbDab 的 1133 条治疗性抗体记录，共 1243 个有效配对条目（含部分双特异性第二臂），不再局限于 8 条 IL4RA 参考。237 个原审查池均完成比较，未新增触发“六 CDR 相似度 ≥0.8 且 H3 相似度 ≥0.8”的候选。
- 相似度为 1 − Levenshtein 距离 / 最大序列长度；六 CDR 按 H1|H2|H3|L1|L2|L3 拼接。全 40 个之间最大六 CDR 相似度 0.784314。
- 参考库采用已冻结本地 CSV（临床状态列标注 Feb 2025），不是实时全库、全部专利或未公开序列。通过本次比对不等于主办方认定原创。

## 完整 IgG 重轻链装配

采用官方模板隐藏表的 Human IgG1（CH1–hinge–CH2–CH3，330 aa）和 Human Ig Kappa（CL，107 aa）。恒定区氨基酸逐字匹配 UniProt P01857-1 的分泌型和 P01834，并独立核对模板中的恒定区 DNA 翻译结果。P01857 的膜型不是本次装配来源。

| 项目 | 本次核验结果 |
|---|---|
| 重链 | VH + 330 aa 恒定区；完整 HC 为 449 或 453 aa |
| 轻链 | VL + 107 aa kappa 恒定区；完整 LC 为 214 aa |
| 配对 | 每候选为 H₂L₂：两条相同 HC、两条相同 LC，无串联 linker、标签或附加信号肽 |
| 接头 | VH/CH、VL/CL 无重复或漏接；逐条接头序列在 JSON |
| Cys | 每条 HC 11 个、LC 5 个；完整 H₂L₂ 的 32 个 Cys 均被预期 16 对二硫键连接表覆盖一次 |
| 糖基化序列 | 每条 HC 保留预期 Fc NST 基序；所选 VH/VL/LC 不含额外 N-X-S/T |
| 模板填写 | C/F 栏只填 VH/VL；D/G 选择 Human IgG1 / Human Ig Kappa；E/H 保留公式与正确缓存 |
| 终止符 | 官方模板恒定区保留原始 **；成熟蛋白 FASTA 删除末尾 **，仅含 20 种标准氨基酸 |
| 导出一致性 | 两份 Excel、80 条 FASTA 链、CSV、JSON 逐条一致；已核验源 CIF SHA256 |

二硫键连接表是基于保守序列的预期拓扑，不代表实验中已经正确氧化成键。位点坐标均是导出完整链的 1-based 连续编号，不是 EU 编号。

## 结构复核范围与限制

对 237 个候选检查 Fv 的 H/L 原子接触与链内 SG 距离，并把 5I19 的配对框架刚性叠合到候选 Fv，检查其已解析的 CH1/CL 常数结构与完整 canonical IL4RA ECD 的冲突。本轮没有发现此固定 Fab 构象的新增恒定区–靶点 <2 Å 冲突。

这只是一个已知 Fab 肘角的局部空间检查。5I19 不含完整 hinge/Fc 与第二臂，且末端链间 Cys 坐标不完整；本轮没有进行完整 H₂L₂ 三维预测、松弛、Fc/第二臂空间采样、显式糖链或膜环境建模。完整 IgG 折叠、表达、二价结合、亲和力及特异性均未被本轮验证。

优先 20 个中 OF3 状态：{'not_evaluated': 15, 'evaluated_uncalibrated': 5}。所有原始 OF3 结果和严格通过标签保留在 JSON，不能把新优先队列解读为原严格门槛通过。

## 优先名单

| 顺序 | 候选 | 生成条件 | HC/LC 长度 | 其他序列提示 |
|---|---|---|---|---|
| 1 | R2_b04_d010 | il4ra_p1_h3_12 | 449/214 | oxidation_susceptible_residue:1 |
| 2 | R1_b01_d012 | il4ra_p1_h3_16 | 453/214 | oxidation_susceptible_residue:1 |
| 3 | R2_b10_d003 | il4ra_p2_h3_12 | 449/214 | oxidation_susceptible_residue:1 |
| 4 | R1_b15_d019 | il4ra_p2_h3_16 | 453/214 | oxidation_susceptible_residue:2 |
| 5 | R1_b12_d001 | il4ra_p1_h3_12 | 449/214 | oxidation_susceptible_residue:2 |
| 6 | R2_b05_d010 | il4ra_p1_h3_16 | 453/214 | oxidation_susceptible_residue:4 |
| 7 | R1_b10_d012 | il4ra_p2_h3_12 | 449/214 | oxidation_susceptible_residue:1 |
| 8 | R2_b11_d013 | il4ra_p2_h3_16 | 453/214 | oxidation_susceptible_residue:1 |
| 9 | R1_b08_d010 | il4ra_p1_h3_12 | 449/214 | oxidation_susceptible_residue:2 |
| 10 | R2_b05_d004 | il4ra_p1_h3_16 | 453/214 | oxidation_susceptible_residue:2 |
| 11 | R2_b03_d006 | il4ra_p2_h3_16 | 453/214 | oxidation_susceptible_residue:3 |
| 12 | R1_b00_d013 | il4ra_p1_h3_12 | 449/214 | oxidation_susceptible_residue:3 |
| 13 | R2_b09_d020 | il4ra_p1_h3_16 | 453/214 | oxidation_susceptible_residue:1 |
| 14 | R1_b15_d016 | il4ra_p2_h3_16 | 453/214 | oxidation_susceptible_residue:1 |
| 15 | R1_b08_d000 | il4ra_p1_h3_12 | 449/214 | oxidation_susceptible_residue:2 |
| 16 | R1_b05_d004 | il4ra_p1_h3_16 | 453/214 | cleavage_motif:1, oxidation_susceptible_residue:1 |
| 17 | R1_b03_d024 | il4ra_p2_h3_16 | 453/214 | oxidation_susceptible_residue:2 |
| 18 | R1_b16_d010 | il4ra_p1_h3_12 | 449/214 | oxidation_susceptible_residue:2 |
| 19 | R2_b17_d005 | il4ra_p1_h3_16 | 453/214 | oxidation_susceptible_residue:3 |
| 20 | R1_b16_d003 | il4ra_p1_h3_12 | 449/214 | oxidation_susceptible_residue:2 |

## 文件与复现

- `exports/priority20_official_template.xlsx`：20 个优先候选，供最终人工查看后用于平台导入；尚未上传。
- `exports/priority20_full_HC_LC.fasta`、`priority20_review.csv`、`priority20.json`：完整成熟链、简表、详细证据。
- `exports/reserve20_*`：20 个备选的同类文件，与优先名单分开，不能一起作为 40 个提交。
- `reports/all1000_reviewed.json`：1000 个设计的原始和新增标签。`policy.json`、`summary.json`、`validation.json`：规则、统计、交叉校验。
- `reports/geometry.json`、`numbering.json`、`numbered_reference_cdrs.json`：局部结构、链编号和参考比较依据。
- `code/number_sequences.py` → `code/review_geometry.py` → `code/finalize_review.py` → `code/validate_exports.py`：使用 RFD4-Proteina `.venv/bin/python` 运行。第一步可从编号缓存续跑，后续覆盖本次复核导出文件，不覆盖两轮生成结果或原筛选报告。

## 来源与方法说明

- [GenScript 官方模板](https://www.genscript.com/gsfiles/techfiles/TurboCHO_IgG_AA_Seperate_Region.xlsx)：2026-11-18 成功下载；原始文件与 SHA256 均保留。
- [UniProt P01857-1](https://rest.uniprot.org/uniprotkb/P01857-1.fasta)、[P01834](https://www.uniprot.org/uniprotkb/P01834/entry)：恒定区独立核对。
- [ANARCII](https://www.nature.com/articles/s42003-026-10186-z)：本轮使用本地已安装版本进行 IMGT 编号。
- [Thera-SAbDab 数据库论文](https://academic.oup.com/nar/article/48/D1/D383/5573951)：本轮使用本地冻结快照作为参考范围，并未宣称穷尽全部已知抗体。
- [LAP 原始研究](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011881)：支持把常见序列基序作为风险线索；本轮为自定义基序扫描，不是运行完整 LAP，也没有套用其全部规则。连续疏水长度 5、CDR 净 K/R-D/E 绝对值 6、两档偏好均为本项目启发式选择，未经过本任务实验校准。
