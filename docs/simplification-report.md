# carta-zarr 全庫簡化報告

日期：2026-09-17 · 範圍：`include/`（1,100 行）、`src/`（8,300 行）、`tests/`（7,900 行）、`CMakeLists.txt`
方法：逐檔通讀 `src/` 與 `include/` 全文，`tests/` 與 CMake 做結構性掃描。

調查當下未動任何程式碼。其後 Tier 1 全部、2.1、2.2、2.4、3.1、3.2、3.4、3.6 已在分支 `sweep-up-after-the-moves` 上套用，各項標有狀態；其餘仍是建議。

---

## 總評

這份程式碼的品質在同類函式庫裡屬於高標。值得先說清楚，因為底下列出的東西容易讓人誤以為相反：

- **模組邊界是想過的，不是長出來的。** `PixelSource`、`Transport`、`SchemaProfile` 三道 seam 各自有 ADR 解釋為什麼在那裡，而且 `src/reduce/` 與 `src/read/` 真的能在沒有 `Store`、沒有 TensorStore 的情況下編譯 — 這不是宣稱，是 CMake 裡那幾個 `INTERNAL` target 在強制執行。
- **註解寫的是「為什麼」而不是「做什麼」。** `chunk_blocks.h` 的 64 MiB 預算、`growing_histogram.h` 的 `[[gnu::noinline]]`、`block_emit.h` 的相對／絕對 channel 編號陷阱，都帶著實測數字。這類註解不能因為「簡化」而刪掉。
- **`pure` 與 `impure` 已經切開。** `pass.cc`／`pass_read.cc`、`pieces.cc`／`pieces_read.cc` 是同一個手法的兩次套用。

所以底下沒有「重寫某個模組」這種建議。真正剩下的是三類：**確定的死碼與筆誤**、**同一個事實寫在多處**、**單一函式內部的樣板噪音**。

熱路徑（`AccumulateRow`、`GrowingHistogram::Add`、`SlabWalk::Over` 的 visitor）**完全不在建議範圍內**。ADR 0005 已經量過那裡的 inlining 值多少，任何把 template 參數換成 `std::function`、或把 per-pixel 迴圈抽成函式的動作都要視為效能變更而非簡化。

---

## Tier 1 — 確定項目（低風險，可直接改）

> **狀態：已全部套用**（2026-09-17，五個 commit）。底下保留原始描述作為紀錄。

### 1.1 兩個 include guard 與檔案路徑不符

| 檔案 | 現有 guard | 應為 |
|---|---|---|
| [src/axis_map.h:8](src/axis_map.h:8) | `CARTA_ZARR_SRC_REDUCE_AXIS_MAP_H_` | `CARTA_ZARR_SRC_AXIS_MAP_H_` |
| [src/store_pixel_source.h:8](src/store_pixel_source.h:8) | `CARTA_ZARR_SRC_REDUCE_STORE_SLAB_SOURCE_H_` | `CARTA_ZARR_SRC_STORE_PIXEL_SOURCE_H_` |

兩個都是檔案從 `src/reduce/` 搬上來時留下的。第二個連檔名都沒跟上（`store_slab_source` → `store_pixel_source`）。全庫只有這兩個不一致，其餘 40 個 header 的 guard 都與路徑相符。

### 1.2 `FilesystemTransport::ListNodes` 裡一段兩條分支完全相同的判斷

[src/zarr/transport.cc:119-126](src/zarr/transport.cc:119)：

```cpp
std::error_code metadata_error;
if (iterator->is_directory() && !std::filesystem::is_regular_file(metadata_path, metadata_error)) {
    if (metadata_error) {
        continue;
    }
    continue;
}
```

`if` 成立與不成立都 `continue`。等價於 `if (...) continue;`，`metadata_error` 也就不需要宣告。這不改變任何行為。

### 1.3 `ProbeReport::ClearDiagnostics` 無人呼叫

[src/schema/xradio/probe_report.h:60](src/schema/xradio/probe_report.h:60) 與 [src/schema/xradio/probe_report.cc:33](src/schema/xradio/probe_report.cc:33)。全庫（含 `tests/`）搜尋只有宣告與定義兩處。`SetDiagnostics` 有一個呼叫者（`InspectImages`），`ClearDiagnostics` 零個。可刪。

### 1.4 `spectral_reduce.cc` 一行孤兒註解

[src/reduce/spectral_reduce.cc:35](src/reduce/spectral_reduce.cc:35)：

```cpp
// Where each axis role sits in the logical axis order.

// One region as the walk sees it: ...
```

它描述的是 `AxisMap`，而 `AxisMap` 已經搬到 [src/axis_map.h](src/axis_map.h)。現在懸在 `WalkRegion` 的註解上面，讀起來像是 `WalkRegion` 的第一句。可刪。

另外 [src/reduce/spectral_reduce.cc:350-355](src/reduce/spectral_reduce.cc:350) 有一段 `// The occupied runs of each chunk row. ... rectangle that happens to contain them.` 的註解，緊接著（[:356](src/reduce/spectral_reduce.cc:356)）卻是 `FootprintBounds` 的註解與定義 — 前一段描述的是更下面的 `BuildColumnRuns`，中間被 `FootprintBounds` 插隊了，兩塊註解之間連空行都沒有。應該把 `FootprintBounds` 移到那段註解之前，或把註解移到 `BuildColumnRuns` 頭上。

### 1.5 `AttributeNumber` 的名字與同組函式不一致

[src/schema/xradio/attributes.h:42](src/schema/xradio/attributes.h:42)：

```cpp
inline bool        HasAttribute   (const json& attributes, std::string_view name);
inline std::string AttributeString(const json& attributes, std::string_view name);
inline std::optional<double> AttributeNumber(const json& value);   // ← 簽章不同
```

前兩個取 `(物件, 名稱)`，第三個取的是**已經取出來的值**。只有一個呼叫者（[image.cc:183](src/schema/xradio/image.cc:183)），叫 `AsNumber` 或 `NumberValue` 會誠實得多。

---

## Tier 2 — 同一個事實寫在多處

這是本庫目前最值得處理的一類。每一項都是「加一個新型別／新欄位時要記得改 N 個地方」。

### 2.1 Zarr data type 表存在四份，各自涵蓋不同子集 ★最高優先

> **狀態：已套用**。新增 [src/zarr/data_type.h](src/zarr/data_type.h)（header-only，不含 JSON 相依），
> 五處改為查同一張表。唯一的行為差異：`ParseDataType` 現在對 complex 回報 `complex64`／`complex128`
> 而非 `unknown`，`DescribeImage` 在此之前就已拒絕非實數型別，故實際不可達。

| 位置 | 形式 | 涵蓋 |
|---|---|---|
| [array_metadata.cc:47 `IsRealDataType`](src/zarr/array_metadata.cc:47) | `array<string_view,11>` | 11 個實數型別，無 bool、無 complex |
| [array_metadata.cc:53 `ParseDataType`](src/zarr/array_metadata.cc:53) | `map<string_view, DataType>` | 12 個，有 bool，無 complex |
| [dataset_size.cc:19 `ElementSizeBytes`](src/dataset_size.cc:19) | `map<string_view, uint64_t>` | 14 個，含 complex |
| [pixel_reader.cc:31 `MatchesDataType`](src/zarr/pixel_reader.cc:31) | 12 條 `if` 鏈 | 12 個，對到 `tensorstore::dtype_v` |

再加上 [chunk_blocks.h:73 `DecodedChunkBytes`](src/chunk_blocks.h:73) 用 `DataType` enum 又寫了一次 bytes-per-element 的 `switch`（含 `default: 4`，意味著 `unknown` 會被當成 4 bytes 而不是被拒絕）。

**建議**：在 `zarr/array_metadata.h` 放一張表，一列一個型別：

```cpp
struct DataTypeInfo {
    std::string_view name;
    DataType         kind;
    std::uint64_t    element_bytes;
    bool             real;          // IsRealDataType 問的那件事
};
```

`IsRealDataType`、`ParseDataType`、`ElementSizeBytes` 三個都變成查這張表。`MatchesDataType` 因為要碰 `tensorstore::dtype_v`（模板，不能進表）留在原地，但可以改成 `switch (ParseDataType(expected))`，於是它至少跟表共用「名稱 → kind」這一半。`DecodedChunkBytes` 改查 `element_bytes`。

**注意**：`chunk_blocks.h` 目前把未知型別當 4 bytes 是刻意的 fallback（讀取預算算錯只是變慢，不是變錯），改成查表時要保留這個 fallback，不要順手改成報錯。

### 2.2 `Image` 五個進入點的樣板完全相同

> **狀態：已套用**。`carta_zarr.cc` 內新增 file-local 的 `WithReadableImage`；`ReadBeams` 維持原樣。

[src/carta_zarr.cc:202, 221, 240, 259, 278](src/carta_zarr.cc:202) — `Read`、`ReadPixelMask`、`ReduceSpectral`、`ComputeHistogram`、`ComputeCubeHistogram` 逐字重複：

```cpp
const std::string node = _impl ? _impl->descriptor.id : std::string{};
return Guarded(ErrorCode::io_error, node, [&]() -> Result<T> {
    if (!_impl /* 或 !_impl || !_impl->store */) {
        return Error{ErrorCode::invalid_argument, "Image handle is empty"};
    }
    auto image = _impl->Readable();
    if (!image) { return image.error(); }
    return internal::X(image.value(), ...);
});
```

而且五份**不一致**：`Read` 與 `ReadPixelMask` 只檢查 `!_impl`，另外三個檢查 `!_impl || !_impl->store`。由於 `Image::Impl` 的 constructor 一定會拿到一個 store，第二個條件永遠為假 — 兩種寫法目前等價，但讀者得自己推導出這件事。

**建議**：加一個私有 helper，把 guard、空 handle 檢查、`Readable()` 一次做完：

```cpp
template <typename Function>
auto Image::WithReadable(Function&& function) const -> decltype(function(std::declval<const internal::ReadableImage&>()));
```

五個進入點各自縮成兩三行，`!_impl->store` 的分歧也就只有一個答案。`ReadBeams` 用不同的 error code（`invalid_metadata`）而且不需要 `ReadableImage`，維持原樣即可。

同樣地 [carta_zarr.cc:187, 192, 362](src/carta_zarr.cc:187) 有三份 `static const T empty; return _impl ? _impl->x : empty;` — 這個重複很小，留著也可以。

### 2.3 三份「把 rows 切成 tasks 丟給 WorkPool」的迴圈

- [plane_histogram.cc:161-183](src/reduce/plane_histogram.cc:161)（每個 plane 切 rows，寫進 `partials` 的一列）
- [plane_histogram.cc:393-408](src/reduce/plane_histogram.cc:393)（跨 plane 切 rows，寫進 `accumulators[task]`）
- [spectral_reduce.cc:735-758](src/reduce/spectral_reduce.cc:735)（切 units，寫進 `partials` 的一段）

三份都是：算 `max_tasks`（pool size ∩ 某個記憶體上限）→ `PlanRowTasks` → `tasks <= 1` 就原地跑 → 否則清零 partials、`workers.Run`、最後把 partials 合回去。

連常數都散在三處，其中兩處還同名同值：

| 位置 | 名稱 | 值 |
|---|---|---|
| [plane_histogram.cc:116](src/reduce/plane_histogram.cc:116) | `kLeastPixelsPerTask` | `1U << 16U` |
| [plane_histogram.cc:278](src/reduce/plane_histogram.cc:278) | `kLeastPixelsPerTask` | `1U << 16U`（同一個檔案裡的第二份） |
| [spectral_reduce.cc:504](src/reduce/spectral_reduce.cc:504) | `kLeastPixelsPerUnit` | `1U << 16U`（同值，不同名） |

另外 `kPartialBudgetBytes` 在 [plane_histogram.cc:114](src/reduce/plane_histogram.cc:114)（64 MiB）與 [spectral_reduce.cc:735](src/reduce/spectral_reduce.cc:735)（16 MiB）同名但**不同值** — 兩邊各有理由，但同名會讓人以為它們該一致。搬到 `tuning.h` 時要順便改名（例如 `kHistogramPartialBudgetBytes` / `kSpectralPartialBudgetBytes`），別把兩個數字併成一個。

**建議**：抽一個 `work_pool.h` 旁邊的 template helper，例如

```cpp
template <typename Body>   // Body(first, last, slot_index)
void OverRowRanges(WorkPool& workers, std::uint64_t row_pixels, std::uint64_t rows,
                   std::size_t max_tasks, std::uint64_t least_pixels, Body&& body);
```

負責 `PlanRowTasks` + 原地／並行的分岔 + 邊界計算，呼叫端只留「清零」與「合併」— 那兩件事三處各不相同（histogram 是加、spectral 有 min/max），不該被一起抽掉。

**必須是 template，不能是 `std::function`** — 理由同 ADR 0005。上表那些常數則提到 `reduce/tuning.h`，那正是 tuning.h 存在的用途。

### 2.4 `nlohmann::json` 的成員探測樣板

> **狀態：已套用**。`ObjectMember` 更名為 `Member`，新增 `MemberObject`／`MemberArray`／`MemberNumber`；
> `src/` 中 `!= nullptr && ->is_X()` 的兩段式寫法歸零。另有四處改用既有的 `AttributeString`。

全庫有 21 處 `ObjectMember(x, "n") != nullptr && ptr->is_object()` 這種兩段式寫法（[image.cc](src/schema/xradio/image.cc) 10 處、[observation.cc](src/schema/xradio/observation.cc) 11 處、[direction.cc](src/schema/xradio/direction.cc) 10 處），另有 16 處 `contains(...) && at(...).is_X()` 的鏈式判斷。

現有的 [`ObjectMember`](src/schema/xradio/attributes.h:47) 只檢查**父物件**是不是 object，不檢查取出來的成員，所以每個呼叫端都得補一次。

**建議**：在 `attributes.h` 補三個回傳 `const json*`（不合條件就 `nullptr`）的取用器：

```cpp
const nlohmann::json* MemberObject(const nlohmann::json&, std::string_view);
const nlohmann::json* MemberArray (const nlohmann::json&, std::string_view);
const nlohmann::json* MemberNumber(const nlohmann::json&, std::string_view);
```

於是

```cpp
if (const auto* t = ObjectMember(image.attributes, "telescope"); t != nullptr && t->is_object()) {
```

變成

```cpp
if (const auto* t = MemberObject(image.attributes, "telescope")) {
```

`observation.cc` 的 `telescope` 區塊（[observation.cc:26-52](src/schema/xradio/observation.cc:26)）收益最大：那個 8 項的 `&&` 條件會縮成 3–4 項，而且不會弱化它刻意逐項檢查、避免 nlohmann 丟例外的那個理由。

同一組取用器也能整理 [`ProbeReport::RequireCoordinateSystem`](src/schema/xradio/probe_report.cc:107) 的五段 `!contains || !is_X || ...` 鏈。

### 2.5 `Store` 裡兩個函式共用同一段前置

[`ReadStringArray1DUncached`](src/store.cc:399) 與 [`ReadStorageLayout`](src/store.cc:424) 開頭都是「讀 node metadata → 讀 array metadata → 兩個都檢查」。可抽成一個回傳 `Result<std::pair<const json&, const ArrayMetadata&>>` 的私有 helper，或更直接地抽成一個只回傳兩個指標的小 struct。收益中等，順手做即可。

### 2.6 `WorkPool` 的取 task 迴圈寫了兩次

[work_pool.cc:68-78](src/work_pool.cc:68)（worker 執行緒）與 [work_pool.cc:119-129](src/work_pool.cc:119)（呼叫端執行緒）是同一個 claim-and-run 迴圈。抽成私有 `void DrainTasks(const Body&, std::size_t worker)` 後兩處各剩一行。這段是並行程式碼，改動要小心，但兩份文字目前完全對稱，抽取是機械的。

---

## Tier 3 — 單一函式內的可讀性

### 3.1 `ReadInto` 的中段驗證區塊

> **狀態：已套用**，抽成 `VerifyStoreMatchesSelection`。

[src/zarr/pixel_reader.cc:111-143](src/zarr/pixel_reader.cc:111)：連續五段「TensorStore 說的 rank / shape / dimension names / data type 要和 canonical metadata 一致」的檢查，夾在「開 array」與「建 index transform」之間，把一個本來就長的函式撐到 150 行。

抽成 `Result<void> VerifyStoreMatchesMetadata(const tensorstore::TensorStore<>&, const PixelSelection&, std::string_view expected_data_type, std::string_view node)`，`ReadInto` 的主線就會變成「檢查請求 → 開 array → 驗證 → 切片 → 轉置 → 轉型 → 讀」，七步各一眼。純搬移，無行為變更。

### 3.2 `BuildChunkBuckets::for_each_cell` 的兩個分支高度重複

> **狀態：已套用**，抽成 `scan_rows`，每個分支只提供 `mark_row`。raster 分支原本「已標記的欄位就跳過像素掃描」
> 的優化保留了下來（改用 `marked` 述詞），`build-release` 的 `pass_timing` 無退化。

[src/reduce/spectral_reduce.cc:120-197](src/reduce/spectral_reduce.cc:120)：`runs != nullptr` 與 `mask != nullptr` 兩條路徑各寫一次外層骨架：

```cpp
const std::uint64_t columns = cx1 - cx0 + 1;
occupied.assign(columns, 0);
for (auto cy = cy0; cy <= cy1; ++cy) {
    std::fill(occupied.begin(), occupied.end(), 0);
    std::uint64_t found = 0;
    const auto y_first = std::max(region.v_start, cy * chunk_v);
    const auto y_last  = std::min(region.v_start + region.v_size, (cy + 1) * chunk_v);
    for (std::uint64_t y = y_first; y < y_last && found < columns; ++y) {
        /* ← 只有這裡不同 */
    }
    flush(cy, occupied);
}
```

把中間那層抽成一個 `mark_columns_of_row(y, occupied, found)` lambda（runs 版與 mask 版各一個），外層骨架就只寫一次，約省 35 行。這段不在 per-pixel 熱路徑上（每個 region 每個 chunk row 走一次），抽取是安全的。

### 3.3 `DescribeImage` 的座標讀取：宣告與使用分離

[src/schema/xradio/image.cc:348-394](src/schema/xradio/image.cc:348)：四個 `std::vector<double>` 在函式中段一次宣告完，然後各自在下面十幾行後才被 `std::move` 賦值。`time_values` 從宣告到使用隔了 45 行，而且只被用一次。

改成在使用處宣告即可：

```cpp
auto l = ReadNumericCoordinate(store, "l");
if (!l) { return l.error(); }
```

四段各省兩行，而且讀者不必記住「上面宣告的那四個 vector 現在填到第幾個了」。

### 3.4 `DescribeSpectralCoordinate` 重複取出同一個成員

> **狀態：已套用**（與 2.4 同一個 commit）。

[src/schema/xradio/image.cc:163-180](src/schema/xradio/image.cc:163)：`ObjectMember(frequency_attributes, "reference_frequency")` 被取了兩次（一次為了 `attrs`，一次為了 `data`），中間只隔了 `spectral.unit` 的賦值。提到區塊開頭取一次即可。

### 3.5 `store_context.cc` 的五個 `NOLINTNEXTLINE`

[src/zarr/store_context.cc:85, 101, 104, 108, 112](src/zarr/store_context.cc:85)：五處都是為了 `spec["key"] = ...` 這個寫法而壓制同一條 clang-tidy 規則。改用 `spec.emplace("cache_pool", ...)` 或 `spec.push_back({...})` 可以讓五個抑制註解一起消失。

順帶一提，`WithoutCache()`（[:84](src/zarr/store_context.cc:84)）與 `MakeStoreContext()` 的 `disable_cache` 分支（[:103](src/zarr/store_context.cc:103)）建的是同一段 `{"cache_pool": {"total_bytes_limit": 0}}`，可共用一個小函式。

### 3.6 `IsNumeric` 是一層沒有內容的包裝

> **狀態：已套用**。

[src/zarr/array_metadata.cc:19](src/zarr/array_metadata.cc:19)：`bool IsNumeric(const json& v) { return v.is_number(); }`，唯一用途是餵給 `std::all_of`。直接寫 `[](const auto& v) { return v.is_number(); }` 即可省掉一個 anonymous namespace。

---

## Tier 4 — 需要先驗證再動的事項

### 4.1 `FilesystemTransport::ListNodes` 的檔案分支可能只產生重複項

[src/zarr/transport.cc:96-103](src/zarr/transport.cc:96) 對「名為 `zarr.json` 的一般檔案」也會 push 它的 parent。但推導下來：

- 對 **array 目錄**：目錄項先 push，接著讀 metadata 發現 `node_type == "array"` 就 `disable_recursion_pending()`，所以它底下的 `zarr.json` **永遠不會被走訪**。
- 對 **group 目錄**：目錄項 push 一次，遞迴進去後它的 `zarr.json` 檔案再 push 一次 → **同一個節點進 list 兩次**。
- 對 **root 的 `zarr.json`**：兩條路徑都被 `relative_parent == "."` 擋掉。

也就是說檔案分支只會產生重複，而 [`Store::ListNodes`](src/store.cc:308) 的 `sort` + `unique` 把它吃掉了，所以外部看不出來。

**若推導成立**，刪掉檔案分支可以讓這個函式少掉一半的條件判斷（`is_directory()` 的三元運算式全部消失）。**但請先驗證**：加一個「有 group 子目錄的 fixture，listing 前不去重」的測試確認重複確實存在，再動。權限錯誤（`skip_permission_denied`）下 `is_directory(error)` 失敗的路徑我沒有完整推完。

### 4.2 `ArrayDirectory` 是 `NormalizeNodeName` 的第二套規則

[transport.cc:170-183](src/zarr/transport.cc:170) 對 node 名稱做的檢查，與 [store.cc:50 `NormalizeNodeName`](src/store.cc:50) 大致重疊但**不完全相同**（`ArrayDirectory` 拒絕 `.`，`NormalizeNodeName` 接受並移除它）。`store.cc` 的註解記載了這個分歧曾經造成 `flag: "./MASK_0"` 的 bug，修法是讓 `Store` 先正規化。現在 `Transport` 的檢查是純防禦性的第二道 — `transport.h` 的 doc 也說「A Transport may check again; it must not have a second opinion」。

目前兩者不會衝突（因為 `Store` 一定先正規化）。但這仍是同一條規則的兩份實作。若要合併，把 `NormalizeNodeName` 移到 `transport.h` 旁邊、讓 `ArrayDirectory` 直接呼叫它是最乾淨的作法 — 這是跨 seam 的改動，值得單獨一個 commit 和一組測試。

### 4.3 `SchemaProfile::Describe` 與 `Discover` 只有測試在用

`Describe` 三個呼叫者、`Discover` 五個，**全部在 [tests/schema_profile_test.cc](tests/schema_profile_test.cc)**。production 路徑走的是 `DescribeVerified`（`Dataset::OpenImage`）與 `entry.inspect`（`ProbeStore`）。

這不一定是問題 — `profile.h` 的註解明說 `Describe` 是「先 inspect 再問能不能開」的版本，是 profile 介面的一部分。但如果要縮小 `SchemaProfile` 的表面積，這兩個是候選。**建議：留著**，只是在文件裡標明它們目前的唯一消費者是測試。

---

## 明確不建議動的地方

為了讓後續的人（或 agent）不要反覆提同樣的建議，這裡列出「看起來像重複，但不該合併」的東西：

1. **`Result<T>` 的 `if (!x) return x.error();`** — 全庫約 200 處。用巨集（`CARTA_ZARR_TRY`）可以砍掉幾百行，但會讓 control flow 藏進巨集裡，對一個給外部 C++ consumer 用的函式庫來說是淨損失。維持現狀。
2. **per-pixel 迴圈的重複**（`AccumulateRow` 的四個 template 實例化、`bin_rows`、`take_rows`）— ADR 0005 量過。不要合併，不要改成 `std::function`。
3. **`RunPass` / `SlabWalk::Over` / `BlockEmitter::Over` 的三層 template** — 已經是上一輪抽取的結果，`block_emit.h` 的註解記載了抽取前的狀態。再往上抽會變成 reduction base class，那正是註解裡說明不要做的事。
4. **`chunk_blocks.h` 與 `tuning.h` 裡帶實測數字的常數註解** — 那些數字是這個庫最難重建的資產。
5. **CMakeLists.txt 的粒度** — 26 個測試 target 看起來冗長，但每一個都有註解說明它為什麼要獨立編譯（多半是為了證明某個模組不依賴 TensorStore）。唯一的小建議：`src/schema/xradio/*.cc` 那組來源在 4 個 target 裡重複列出，可以比照現有的 `CARTA_ZARR_METADATA_SOURCES` 再設一個變數。

---

## 測試端

`tests/` 沒有逐檔通讀，但結構掃描出一項明確的重複：

**24 個測試檔各自定義一份完全相同的 `Require`**：

```cpp
void Require(bool condition, const std::string& message) {
    if (!condition) { throw std::runtime_error(message); }
}
```

此外 `RequireClose` 有兩份不同簽章（[conformance_test.cc:38](tests/conformance_test.cc:38) 帶 tolerance 參數，[spectral_reduce_test.cc:56](tests/spectral_reduce_test.cc:56) 不帶），而 `main()` 的 harness 也是 24 份手寫，且分成兩種風格 — 有的用 `std::cout` 並在成功時印訊息，有的用 `std::fprintf` 且成功時不印。

**建議**：新增 `tests/support/check.h`（header-only，與現有的 `in_memory_transport.h`、`synthetic_pixel_source.h` 同一個位置），放 `Require`、`RequireClose`，以及一個接受 `{name, fn}` 清單的小 runner：

```cpp
int RunTests(std::string_view suite, std::initializer_list<std::pair<std::string_view, void(*)()>>);
```

省下約 150 行，並且讓失敗輸出格式一致。**注意**：這會讓每個測試 target 的 SOURCES 多一個 header 依賴（header-only 的話不必改 CMake），且要確認那些刻意不連結任何東西的 target（如 `carta_zarr_linear_axis_tests`）不會因此被迫連結新東西。

---

## 建議的執行順序

| 順位 | 項目 | 風險 | 行數變化 |
|---|---|---|---|
| ~~1~~ | ~~Tier 1 全部（1.1–1.5）~~ **已完成** | 無 | −5 淨（−24/+19，含兩處註解搬移與補充） |
| ~~2~~ | ~~2.1 data type 表合一~~ **已完成** | 低 | 消費端 −9；新表 +81（含註解），換來「加型別只改一處」 |
| ~~3~~ | ~~2.4 JSON 取用器 + 3.4、3.6~~ **已完成** | 低 | −8 淨（消費端 −60，取用器 +32） |
| ~~4~~ | ~~2.2 `Image` 進入點 helper~~ **已完成** | 低 | −12 |
| ~~5~~ | ~~3.1、3.2 兩處函式內抽取~~ **已完成** | 低 | 實際 **+38**（新增 15 行註解說明抽取理由；可讀性為主，行數本非目標） |
| 6 | 測試端 `check.h` | 低 | −150 |
| 7 | 2.3 `OverRowRanges` | **中**（碰並行與熱路徑邊緣） | −50 |
| 8 | 2.6 `WorkPool::DrainTasks` | **中**（並行） | −15 |
| 9 | 4.1 transport listing | **中**（需先加測試） | −25 |

1–6 合計約減 200 行，且沒有一項會改變任何可觀察行為。7–9 建議各自獨立 commit，並在 `build-release` 上跑一次 `pass_timing` 對照。

---

## 後續驗證

每一批改動後：

```bash
cmake --build build > /tmp/carta_zarr_build.log 2>&1 || (grep -in -C 3 "error:" /tmp/carta_zarr_build.log | head -n 50 && exit 1)
```

```bash
ctest --test-dir build --output-on-failure > /tmp/carta_zarr_test.log 2>&1 || (tail -n 60 /tmp/carta_zarr_test.log && exit 1)
```

Tier 2 第 7 項與之後的改動另需在 `build-release` 跑 `carta_zarr_pass_timing` 對照。
