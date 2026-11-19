# GenScript Challenge — IL4RA paired IgG computational campaign

两轮各 500 个 de novo 配对抗体设计的计算记录、筛选证据与 IgG 序列交付材料。目标为人 IL4RA（UniProt P24394-1，26–232 aa 胞外域）。两轮计算在集群节点 p-gc-3039 完成；本仓库归档日期为 2026-11-19。

**当前状态：1000 个设计生成完成，最终优先 20 个＋备选 20 个已完成完整重轻链序列装配。尚未提交比赛平台，尚无实验结合结果，也没有完成完整 H₂L₂ 三维预测。原 OF3 严格全门槛通过数仍为 0。**

## 交付文件

- [优先 20 个：已填写的官方 Excel](computational/igg_review_20261118/exports/priority20_official_template.xlsx)
- [优先 20 个：完整 HC/LC FASTA](computational/igg_review_20261118/exports/priority20_full_HC_LC.fasta)
- [优先 20 个：序列与指标 CSV](computational/igg_review_20261118/exports/priority20_review.csv)
- [完整筛选与 IgG 装配复核报告](computational/igg_review_20261118/README.md)
- [结果交付 ZIP](computational/igg_review_20261118/IgG_review_package.zip)
- [计算记录索引](computational/igg_review_20261118/计算记录索引.md)

平台导入应使用 **priority20**，不要把备选合并成 40 个提交。官方表仅在 VH/VL 栏填写可变区，恒定区选择 Human IgG1 / Human Ig Kappa。完整 HC 为 449 或 453 aa，LC 为 214 aa；FASTA 为成熟蛋白序列，不附加信号肽。上传文件后仍需平台校验与订单确认。

## 计算结果概览

| 项目 | 第一轮 | 第二轮 | 合并 / 最终 |
|---|---:|---:|---:|
| 新生成设计 | 500 | 500 | 1000 |
| 原始几何/序列初筛通过 | 219 | 186 | 405 |
| 独立 VH/VL 序列对 | 500 | 500 | 998（跨轮存在重复） |
| 合并去重及额外序列筛选后的审查池 | — | — | 237 |
| 追加序列与配对结构复核后的池 | — | — | 211 |
| 最终队列 | — | — | 优先 20，备选 20 |

优先 20 个均没有本次重点检查的 CDR N[GS]/D[GS] 基序。备选中 9 个同样无这些基序，11 个带标记保留。没有放宽原有几何筛选门槛；仅在补齐备选时放宽化学基序偏好。基序缺失不等于不存在表达、降解或聚集风险。

原始 OF3 评估与诊断续算均保留。已知阳性对照未通过校准，因此后续选择不把 OF3 结果作为可靠否决或加分项；`not_evaluated` 与 `evaluated_uncalibrated` 分开记录。优先队列不能解释为实验结合预测成功，或原严格门槛通过。

## 方法与参数

- RFD4-Proteina revision `aaccdd2ad0720c180468df5f9cf037a713d403e5`；模型快照 `a06e2509c7033b010fc43614e9fd840303588502`。模型权重、外部代码仓库和运行环境不包含在本仓库。
- 采用 5I19 的通用人源 IGHV3-23 / IGKV1-39 配对框架，六条 IMGT CDR 均作为生成区域；已知 IL4RA 抗体仅用于参考比较或诊断对照，不作为生成起点。
- 两轮各 20 批，每批 25；每轮四条件各 125：`p1/p2 × H3 length 12/16`。p1 热点为 107/108/109，p2 为 67/70。
- 第一轮种子 `2026111700–2026111719`；第二轮 `2026111900–2026111919`。生成设置保持一致，仅 seed / output directory 改变：medium、front_loaded、60 diffusion steps、180 trunk passes、BF16、batch size 1。
- 详细阈值与选择顺序见 [最终筛选规则](computational/igg_review_20261118/reports/policy.json)、[原合并规则](computational/combined1000_20261118/policy.json) 与逐批配置。相似性阈值和软风险规则均为项目内部筛选条件，不是主办方官方规则。
- 扩展参考比较覆盖本地冻结 TheraSAbDab 的 1133 条记录 / 1243 个配对条目，非全面专利或全部已知抗体检索。
- 完整链装配使用官方模板恒定区，并与 UniProt 和模板恒定区 DNA 翻译交叉核验。结构复核仅包括 Fv 和固定 Fab 恒定区空间兼容性，不覆盖完整 Fc、第二臂、显式糖链或膜环境。

## 仓库结构

| 路径 | 内容 |
|---|---|
| `Challenge_任务梳理与背景知识.md`、根目录 PDF | 原任务依据与解读 |
| `computational/scripts/`、`configs/`、`provenance/` | 计算脚本、参数与来源记录 |
| `computational/campaign_local500_20261117/` | 第一轮批次参数、运行日志、逐条指标、初次 OF3 评估 |
| `computational/campaign_local500_round2_20261118/` | 第二轮批次参数、解析配置、运行日志与逐条指标 |
| `computational/continuation_20261118/` | 阳性对照、四个候选的复测、失败尝试及修复验证 |
| `computational/combined1000_20261118/` | 初版统一筛选，保留历史名单 |
| `computational/igg_review_20261118/` | 追加复核、最终名单、完整链、Excel 与验证记录 |
| `computational/campaign_36h*/`、`runs/`、`reports/` | 较早试运行、失败记录与准备阶段材料 |
| `archives/` | 去重压缩的原始 CIF、SQLite 台账与大编号 JSON，附逐文件映射与 SHA256 |
| `repository_tools/` | 无第三方依赖的归档、校验与恢复脚本 |

第一轮状态 `completed_with_stop_condition` 是 OF3 停止条件触发，**不表示 500 个生成未完成**。条件式 32-variant 分支没有启动；保留计划文件不代表其已经执行。旧文档中描述的运行中状态属于历史记录，以本页和最终复核报告为准。

## 恢复原始计算文件

原始目录中重复保存的 CIF 按 SHA256 去重，分片压缩存入 Git；没有使用 Git LFS，也不需要额外下载凭据。日志、配置和结果表直接保留为可阅读文件。原始文件路径、哈希和所属归档见 [file_map.csv](archives/file_map.csv)。

```bash
git clone https://github.com/SaaaaiDing/Genscript_Challenge.git
cd Genscript_Challenge
python3 repository_tools/restore_artifacts.py --verify-only

# 恢复到当前 clone；已有文件必须哈希一致，否则拒绝覆盖。
python3 repository_tools/restore_artifacts.py --destination "$PWD"

# 也可只恢复第一轮：
python3 repository_tools/restore_artifacts.py --destination "$PWD" \
  --prefix computational/campaign_local500_20261117/
```

完整恢复会重新展开重复存放的原始文件，需要数 GB 空间。模型权重、环境缓存、临时文件、凭据 `.env`、锁文件和重复的非 CIF 内容寻址副本不上传。源配置与报告副本保留；归档覆盖范围见 [manifest.json](archives/manifest.json)。结果交付 ZIP 仅是交付包，完整计算归档由本 Git 仓库的日志/配置及 `archives/` 共同构成。

## 可复现范围

归档可在仅有 Python 标准库的环境下校验与恢复。重新分析或生成需要原 RFD4-Proteina 环境、ANARCII / Biotite 等依赖、对应模型权重、GPU 与 OF3 运行环境。原脚本、SQLite 和配置中的 `/storage/...` 绝对路径、Slurm job ID、节点名刻意保留为原始执行证据；迁移运行时应先映射这些路径和资源设置，不能在普通 clone 上直接启动原集群脚本。

未附加开源许可证：此归档并不重新授权第三方模型、参考数据、模板或比赛文件。比赛实际提交、实验结果及主办方新颖性判定不由本仓库替代。
