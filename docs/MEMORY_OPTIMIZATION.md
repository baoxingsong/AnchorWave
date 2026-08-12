# AnchorWave 非比对工作区内存优化与评估

## 目标与边界

本轮优化针对序列比对工作区以外、随染色体长度或锚点/SAM 记录数增长的常驻与瞬时内存。WFA、BiWFA 和 KSW2 自身的工作区、`-w` 的含义，以及全进程内存调度策略不在本轮修改范围内。

实现顺序与预定优先级一致：

1. 取消整条染色体比对字符串的内存缓存。
2. 消除锚点和 SAM 容器的重复副本。
3. 压缩 `Transcript` 和 `AlignmentMatch`。
4. 限制 `readSam()` 的评分缓存。
5. 减少 FASTA 区间读取的瞬时副本。

## 实现

### 1. 不再在内存中累计整条染色体比对

`genomeAlignmentSingleThread()` 和 `genomeAlignmentAndVariantCallingSingleThread()` 不再用两个 `std::stringstream` 累积 reference/query 的整条比对结果。

- 不输出 whole MAF 时，只保留长度计数，并逐区间验证去 gap 后的序列，不保存历史比对行。
- 输出 whole MAF 时，reference/query 行分别写入匿名临时文件，最终以 64 KiB 缓冲顺序流式导出。
- 临时文件由 `std::tmpfile()` 管理，关闭时自动清理；所需临时磁盘空间约等于两条带 gap 比对行的总长度。
- 去 gap 一致性检查改为直接扫描，不再构造临时 ungapped 字符串。

因此，whole MAF 关闭时，这部分 RAM 从 `O(染色体累计比对长度)` 降为 `O(当前区间长度)`；whole MAF 开启时，RAM 同样有界，但会用临时磁盘换取内存。

### 2. 去除容器副本

- 染色体级 `AlignmentMatch` vector 进入任务时使用 move，不再逐任务复制。
- SAM 原始匹配结果转入染色体 map 后立即释放原容器。
- 任务之间共享只读 FASTA 索引/序列 map，不再为每个线程复制。
- `Transcript` 的 name index 保存 `const Transcript*`，不再保存完整 `Transcript` 副本。
- 输出、排序和遍历路径尽量使用 const reference；可转移的中间结果使用 move。
- GFF 解析完成后把 `Transcript` 从 name map 移入 chromosome vector，避免最后一次完整复制。

### 3. 压缩 `Transcript` 与 `AlignmentMatch`

- `Transcript` 删除了与 CDS/exon vector 重复、且未被使用的两个 set 成员。
- `AlignmentMatch` 中四个高重复字符串（reference/query 染色体名和基因名）改为进程生命周期字符串池中的指针。
- `AlignmentMatch` 复制变为固定大小字段复制，不再复制四个 `std::string` 的对象与潜在堆内存。
- 字符串池只保留唯一值且不删除，因此适合染色体名和基因名这类高度重复、唯一值规模可控的数据。

单元测试要求 `sizeof(AlignmentMatch) < 96`，并验证相同字符串实际共享地址。

### 4. 限制 `readSam()` 评分缓存

原实现为每个 gene/chromosome 组合保留全部命中评分。当前实现只保留最高的 `expectedCopies + 1` 个评分：前 `expectedCopies` 个代表允许的拷贝，额外一个只用于判断是否存在相似度过高的多余拷贝。

空间复杂度从 `O(全部 SAM 命中数)` 降为 `O(gene/chromosome 组合数 × (expectedCopies + 1))`，判定语义不变。

### 5. 减少序列读取瞬时副本

未压缩 FASTA 的区间读取现在把 `pread()` 结果直接写入已定长的 `std::string`，随后原地移除换行并缩短到准确碱基数。不再经历 `char[] -> string -> substr` 的三阶段复制。FASTA index 参数也改为 const reference。

压缩输入的既有自动识别与读取路径保持不变，并由原有 compressed-input 测试覆盖。

## 测试

构建目录：`/private/tmp/anchorwave-memory-after`

`ctest --output-on-failure` 共 7 项全部通过：

- trio unit
- WFA2 BiWFA unit
- anchor scheduler unit
- compressed input unit
- memory unit
- trio CLI smoke
- trio workflow

新增 memory unit 覆盖：

- alignment spool 的流式写出、长度统计、无 gap 验证和异常输入。
- `AlignmentMatch` 字符串共享与对象尺寸。
- `readSam()` top-K 评分缓存的容量和判定边界。

## B73/Mo17 A/B 评估

> 注：本节 A/B 评估早于 v1.3.0 对 `-fa3^2` 32 位溢出的修复。它用于衡量内存优化和同版本输出一致性；其中 109,413 行及其哈希不是 v1.3.0 修复后的发布金标准。

对照二进制是在本轮五项内存修改前保存的 `/private/tmp/anchorwave-memory-before/anchorwave`；测试二进制是 `/private/tmp/anchorwave-memory-after/anchorwave`。两者使用同一数据、参数和算法选择逻辑。

| 场景 | 线程 | 修改前峰值 RSS | 修改后峰值 RSS | 变化 | 修改前墙钟 | 修改后墙钟 |
|---|---:|---:|---:|---:|---:|---:|
| 全量 B73/Mo17，`genoAli` 仅生成 anchor | 10 | 1,625,360 KiB | 879,920 KiB | **-45.9%** | 153.27 s | 154.57 s |
| chr1 前 4 Mb，`genoAli` 完整序列输出 | 6 | 14,542,832 KiB | 14,447,920 KiB | -0.65% | 53.41 s | 53.85 s |
| chr1 前 4 Mb，`proali` 完整序列输出 | 6 | 19,605,904 KiB | 22,598,912 KiB | +15.3% | 163.76 s | 164.07 s |

输出一致性：

- 全量 anchor：109,413 行；忽略包含输出路径的首行后，SHA-256 均为 `24fa6b9fac565f5860a0379ea3b51f065c9f0388724b5835fd92dc6ad6866ade`。
- 4 Mb `genoAli`：whole MAF、fragment MAF、methods BED 均逐字节一致；anchor 在规范化首行后相同。
- 4 Mb `proali`：whole MAF、fragment MAF、methods BED 均逐字节一致；anchor 在规范化首行后相同。

## 结果解释

全量 anchor-only 测试隔离了本轮目标内存，峰值减少约 728 MiB（45.9%），而运行时间基本不变。这证明消除 `Transcript`、锚点和 SAM 容器副本对真实全基因组输入有效。

4 Mb 序列比对的总峰值主要由多个并发 WFA/BiWFA 工作区在何时重叠决定；固定结构内存只占其中很小一部分。因此 `genoAli` 混合峰值仅小幅下降，而一次 `proali` A/B 运行出现了 15.3% 的负向峰值波动，尽管算法选择统计、输出和耗时均相同。该结果不能证明固定内存增加了约 3 GiB：新引入的 spool 在本例只有约 18 MiB，量级不符；更合理的解释是并发工作区峰值重叠的运行间波动。若要定量评估序列阶段，应分别进行多次重复，或增加“非比对基础 RSS”和“工作区已预留/实际使用内存”的阶段采样。

## 使用注意事项与后续建议

- whole MAF 开启时，应确保临时目录有足够空间；不需要 whole MAF 时不要传 `-o`，此时不会建立 spool。
- 进程级字符串池用互斥锁保护且生命周期与进程相同。若将来出现数千万个几乎不重复的 gene ID，应改为批次级 arena 或整数 ID 表。
- 下一步最有价值的测量改进，是在调度器中记录 `baseline RSS`、每个 alignment task 的开始/结束 RSS，以及每种算法的实际峰值工作区；这样可以把“约 20 GB 固定部分”和“线程数 × 工作区”可靠拆开。
