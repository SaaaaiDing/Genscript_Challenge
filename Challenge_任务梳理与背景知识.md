整理日期：2026-11-17

依据：本目录《2026_Genscript x Geodesic Challenge Brief_Round 1 1.pdf》，共 10 页，封面日期 2026-11-16。比赛要求以下述 PDF 为准；外部资料只用于解释背景，不代表主办方补充规则。本文是阅读整理，尚未开展候选设计或提交。

**任务的核心**

设计最多 20 个全新的、靶向人 IL-4Rα 胞外域的 IgG 抗体，提交每个候选配对的重链和轻链氨基酸序列。主办方进行实验处理和结合筛选，第一轮筛出能结合靶点的候选，第二轮进一步评价可开发性。

**比赛明确要求（PDF 第 2–6 页）**

| 项目 | 要求 |
|---|---|
| 靶点 | 人 Interleukin-4 Receptor Alpha，IL-4Rα / CD124 |
| 序列基准 | UniProt P24394-1，氨基酸 26–232，即 207 个残基的胞外域 ECD |
| 设计对象 | IgG 格式抗体，须有配对重链和轻链 |
| 数量上限 | 每队最多 20 个独特候选；操作示例中的 20 不等于必须交满 20 |
| 方法 | 使用 AI 驱动的设计方法；允许其他计算建模和抗体工程方法，不限定只用 Geodesic |
| 新颖性 | de novo；不能从已有结合抗体修改而来，尤其明确禁止以已有 IL-4Rα binder 为起点改序列 |
| 可用既有数据 | 公共蛋白数据及已知 binder 可用于训练、基准测试和设计评价 |
| 第一轮实验 | Carterra LSA 单浓度（1-con）结合筛选 |
| 第一轮评价 | 可检测结合信号，以及候选间的相对结合响应；阈值和晋级标准另行通知 |
| 截止时间 | 2026-11-23 23:59 EST（纽约时间），即北京时间 2026-11-24 12:59 |
| 结果时间 | 预计 2027-01-13 至 01-18；属于预期时间，不是保证 |
| 第二轮 | 入选候选接受更多生物物理和可开发性表征；本文件没有列出具体完整检测项目 |
| 可能补交 | 方法说明、预测结果、结合界面或结构模型，主办方可能邮件索取 |

新颖性主要应在 CDR 层面比较。PDF 明确承认框架区和恒定区可以是标准化的，因此不应把“全新”理解为整条 IgG 每个位置都必须创新，也不能只改恒定区来制造表面上的序列差异。文件没有给出相似度淘汰阈值。其“冻结参考库”“允许申诉”等文字是实施建议，不能当成已落实的评审程序。（第 3 页）

**先理解这些背景概念**

1. IL-4Rα 是受体，IL-4 和 IL-13 是信号分子。

   IL-4Rα 位于细胞膜上，其胞外部分参与接收信号。IL-4 可通过 IL-4Rα 与共同 γ 链组成的 I 型受体，或 IL-4Rα 与 IL-13Rα1 组成的 II 型受体传递信号；IL-13 先结合 IL-13Rα1，再招募 IL-4Rα，形成 II 型信号复合物。因此 IL-4Rα 是两条信号路线的共享组分，但不能简单理解为 IL-4 和 IL-13 都以相同方式直接结合孤立的 IL-4Rα。[UniProt 官方条目](https://www.uniprot.org/uniprotkb/P24394/entry)

   PDF 将其定位为与特应性皮炎、哮喘等 2 型炎症相关的已验证治疗靶点。本轮的明确任务是“结合”；能结合受体、能阻断信号、最终能成为药物，是不同层面的结论，不能互相替代。

2. IgG 抗体是由重链和轻链装配起来的蛋白。

   典型 IgG 有两条相同重链和两条相同轻链。一个候选只需定义一条重链序列和一条轻链序列，表达时装配为完整分子。

   - 重链：VH（可变区）— CH1 — hinge（铰链）— CH2 — CH3。
   - 轻链：VL（可变区）— CL（恒定区）。
   - VH 与 VL 配对形成识别抗原的部位；两者各有 3 个 CDR（互补决定区），共 6 个 CDR 环。
   - FR（framework，框架区）支撑可变区结构；CDR 通常承担主要识别作用，但实际接触不一定严格限于 CDR。
   - Fc 是 IgG 的恒定区部分，与效应功能等有关；它不等于整个重链恒定区。

   因此，我们主要要解决的是能正确配对、折叠并识别目标的 VH/VL 设计，再按标准构建体配上恒定区。仅提交 CDR-H3、短肽、单独 VH 或任意蛋白 binder，不能满足文件里的配对 IgG 要求。[IMGT 官方 FAQ](https://www.imgt.org/FAQ/)、[IMGT 结构与命名说明](https://pmc.ncbi.nlm.nih.gov/articles/PMC7555362/)

3. Epitope 与 paratope 分别是接触的两面。

   Epitope（表位）是靶蛋白上被抗体识别的区域；paratope 是抗体上识别靶点的区域。结构设计要考虑两面的形状、电荷和相互作用是否合理，以及抗体是否真的能靠近该位置。得到一个看似接触紧密的模型，并不等于实验一定能结合。

   结构驱动的 de novo 抗体设计可以在框架条件下生成新的 CDR，并用复合物建模和实验验证逐步筛选；这是技术路线示例，不是比赛指定算法。[原始研究：结构驱动的表位特异抗体设计](https://www.nature.com/articles/s41467-025-67361-9)

4. LSA、SPR、亲和力与多价效应。

   Carterra LSA 使用表面等离子体共振（SPR）技术，可以并行监测结合响应。实验通常记录分子结合与解离的时间曲线；平台本身能开展动力学和亲和力测量，但本比赛指定的是单浓度筛选。[Carterra LSA 官方介绍](https://carterra-bio.com/lsa/)

   在适用的简单 1:1 模型中，平衡解离常数 K_D = k_off / k_on，K_D 越低通常表示亲和力越强。但响应还受表面固定量、分子活性、计量关系及实验布局等因素影响，不能把单个响应值直接换算为可靠的 K_D。IgG 有两个结合臂，多位点接触还可能产生 avidity（多价结合带来的整体结合增强）。据此解读本题：第一轮主要确认是否有结合并进行相对筛选，不能直接把其名次等同于严格的单价亲和力名次。[Carterra 实验与数据解释 FAQ](https://carterra-bio.com/customer-faqs/)

5. Developability（可开发性）是能否顺利制备并保持合适性质。

   需要考虑表达、稳定性、聚集、自相互作用、非特异结合及化学降解等。计算上可检查异常疏水或电荷表面、潜在序列风险和结构异常。这些是风险筛查，不是实验结论，也没有一种分数能保证全部性质合格。比赛第一轮虽以结合为主，严重的表达或稳定性问题仍可能影响实际检测。[原始研究：Therapeutic Antibody Profiler](https://pmc.ncbi.nlm.nih.gov/articles/PMC6410772/)

**这个靶点建模时要留意什么**

- PDF 引用的 [PDB 3BPN](https://www.rcsb.org/structure/3BPN) 是 IL-4 / IL-4Rα / IL-13Rα1 三元复合物，不是抗体—IL-4Rα 复合物。不能把整个结构中的所有蛋白都当作本轮靶标；应先识别目标链，并与 P24394-1 的 26–232 区间核对。条目记录有工程突变，不能默认其序列等同于比赛抗原。
- 第 7 页要求选择 Fc-tag 的 IL-4Rα 抗原。与名称和区间相符的 ACROBiosystems ILR-H5253 产品为 HEK293 表达，包含 26–232 区间，C 端连接人 IgG1 Fc，且有糖基化。PDF 未明确给货号，因此这是匹配资料，不是对最终实验批次的确认。[厂家产品资料](https://www.acrobiosystems.com/products/protein/ilr-h5253)
- 由上述信息推得：需要考虑糖链是否遮挡候选表位，以及 Fc 标签和多价呈现是否影响结合构象或测量。裸蛋白模型可能没有呈现这些因素；不要把对标签的相互作用误认为靶点特异结合。这是建模与结果解释建议，不是主办方已公布的扣分规则。

**建议执行顺序与应保留的产物**

下表是工作建议，不是官方额外交付要求。当前日期距截止约 6 天，宜先完成账号和输入规范确认。

| 顺序 | 要做的事 | 建议保留的产物 |
|---|---|---|
| 1 | 核对比赛专用 GenScript 账号；使用报名邮箱注册 Geodesic，检查 token；打开模板 | 账号可用性记录、模板原件 |
| 2 | 核对靶序列、ECD 区间、结构链、缺失残基、突变及糖基化信息 | 靶点 FASTA、结构、编号映射与注释 |
| 3 | 选择多个合理表位或设计假设，生成足够多的 de novo 配对 VH/VL 候选 | 候选表、模型版本、参数与生成记录 |
| 4 | 进行序列与结构质控，检查配对、构象、碰撞和表达/稳定性风险 | 每个候选的筛选记录 |
| 5 | 对 CDR 做统一编号、已知抗体相似性比较和候选间去重 | 新颖性及多样性记录；不虚设官方阈值 |
| 6 | 综合筛选最多 20 个候选，保留不同序列家族和合理界面方案 | 最终名单、选择理由、备选名单 |
| 7 | 装配完整重轻链；按模板填写；校验序列；提前提交 | Excel、完整序列备份、方法摘要、订单号 |

可以按 11 月 17 日处理账号与靶点资料、18–20 日生成与初筛、21–22 日复核与整理、23 日截止前完成提交来安排。实际生成规模取决于 token 和算力；PDF 没有要求必须生成多少中间候选。最终不要只按单个模型分数挑选 20 个几乎相同的序列，文件明确鼓励多样性。

**提交格式和完整操作链（PDF 第 6–9 页）**

第 6 页说最终序列需包含可变区与恒定区；第 8 页则要求在分区模板填写 Antibody name、VH Seq、CH Type、VL Seq、CL Type，CH Seq 根据 CH Type 自动填入。合理的操作理解是：最终表达完整 IgG，上传表采用分区输入。应同时保留完整重轻链备份，但不要把整条重链直接填进 VH 栏，避免恒定区重复。

- [官方 Excel 模板](https://www.genscript.com/gsfiles/techfiles/TurboCHO_IgG_AA_Seperate_Region.xlsx)：首次程序访问返回 HTTP 403；2026-11-18 的 IgG 装配复核已成功下载并核验工作簿。原件保存在 `computational/igg_review_20261118/raw/official_template.xlsx`，隐藏表的 Human IgG1 / Human Ig Kappa 恒定区分别为 330 / 107 aa，与 UniProt P01857-1 / P01834 一致；模板末端 `**` 在成熟蛋白 FASTA 中去除，在官方模板恒定区中保留。
- CH Type：PDF 第 8 页允许 Human IgG1 或 Human IgG4；第 7 页示例提示词用 IgG1。不是只允许 IgG1。
- CL Type：应与 VL 的轻链类型匹配。已核验模板包含 `Human Ig Kappa`、`Human Ig Lambda`；当前 IGKV1-39 框架采用 `Human Ig Kappa`。
- 输入氨基酸序列，不是 DNA。不要自行加入信号肽，平台会添加 `MGWSCIILFLVATATGVHS`。
- [提交入口](https://www.genscript.com/ai-drug-discovery-solutions.html?src=mostpopular)：右下角 AIDD Agent / Dr. G，展开后使用比赛专用账号登录并新建会话。

完整流程：

1. 按 PDF 示例提出 IgG 表达及 Human IL-4Rα 的 LSA 单浓度结合检测请求，数量应与实际候选数一致。
2. 选择 `Human IL-4 R alpha / CD124 Protein, Fc Tag, low endotoxin (MALS verified)`，规格 `100 μg/EA`。
3. 核对 Proposal Summary：IgG 格式、实际序列数量、Carterra LSA 检测。示例单序列费用 152.5 USD 是界面标价；比赛用 credits 抵扣。
4. 根据原文回复 `No quote needed, create an order directly.`，随后上传填写好的 Excel。
5. 通过 Data Validation 与后续生信检查；任何 Failed 序列须修正重传，不应继续确认订单。
6. 按 PDF 指定设置 PO `BB20261112`、Project Name `AIDD-test`、No Comments，并要求先展示订单预览。
7. 核对候选、抗原、服务及 Beyond Binding Challenge Credit；最终预计总额必须是 `0.00 USD`。PDF 明确要求非零不要继续。
8. 预览无误后确认创建订单；必须看到 Order Summary，保留 Order No.。仅上传 Excel 不等于提交成功。

PDF 的“20 records”“22 items”对应提交 20 个候选的示例，其中另外两项是抗原和检测服务。少于 20 个候选时，不应把示例条目数当成硬性要求；界面处理方式应据实际记录核对。

**仍未明确的事项**

- LSA 的具体浓度、固定/捕获方式、参考扣除方法、阳性阈值和晋级规则。
- 新颖性参考库与具体相似度判定阈值。
- 完整 H₂L₂ 的实际折叠、表达与结合仍需验证；模板恒定区版本、分区填写和完整序列拼接已于 2026-11-18 复核，见 `computational/igg_review_20261118/README.md`。
- 第二轮完整检测项目和评分权重。
- 最终实验抗原的货号、批次及实际构建体信息。

这些不妨碍先准备资料和设计流程，但不能把未公布内容写成既定规则。比赛规则与实验问题联系人：`proteinchallenge@genscript.com`；Geodesic 账号/token/平台问题联系人：`contact@geodesiclab.com`。PDF 建议邮件注明姓名、联系邮箱，主题含 Beyond Binding Challenge。本次没有发送邮件。

**后续阅读入口**

- [Geodesic 平台](https://www.geodesiclab.com/)、[帮助中心](https://agent.geodesiclab.com/help)：来自 PDF；本次未验证登录后的模型、配额或使用能力。
- [比赛 Terms and Conditions](https://www.genscript.com/genscript-lsrg-x-geodesic-beyond-binding-challenge-terms-and-conditions.html)：PDF 明示提交即同意；本次未进行条款逐项审阅。
