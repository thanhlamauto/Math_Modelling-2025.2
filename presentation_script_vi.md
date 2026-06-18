# Script thuyết trình 15 phút

Bản này bám theo `slide.tex` sau khi `git pull` từ GitHub. Tổng thời lượng được chia đúng 15:00 cho 16 slide.

Lưu ý khi tập:

- Nói chậm vừa phải, không cố đọc quá nhanh ở các slide công thức.
- Nếu bị chậm hơn mốc, bỏ bớt câu giải thích phụ, nhưng không bỏ slide kết luận.
- Đến mốc 14:30 phải chuyển sang slide Conclusion.
- Nếu chia 5 người: mỗi người nhận khoảng 3 slide, người cuối nhận 4 slide ngắn.

## Bảng phân bổ thời gian

| Slide | Nội dung | Thời gian | Mốc kết thúc |
|---|---|---:|---:|
| 1 | Title | 0:30 | 0:30 |
| 2 | Roadmap | 0:40 | 1:10 |
| 3 | Problem | 1:00 | 2:10 |
| 4 | Modeling Approach | 1:10 | 3:20 |
| 5 | Two Superspreader Models | 1:20 | 4:40 |
| 6 | Infection Kernels | 1:10 | 5:50 |
| 7 | R0 and Critical Density | 1:10 | 7:00 |
| 8 | Monte Carlo Algorithm | 1:00 | 8:00 |
| 9 | Result 1: Percolation | 1:15 | 9:15 |
| 10 | Result 2: Propagation Speed | 1:00 | 10:15 |
| 11 | Result 3: Epidemic Curves | 0:55 | 11:10 |
| 12 | Result 4: Secondary Infections | 0:55 | 12:05 |
| 13 | Comparison with SARS Data | 0:55 | 13:00 |
| 14 | Sensitivity and Robustness | 0:55 | 13:55 |
| 15 | Limitations | 0:35 | 14:30 |
| 16 | Conclusion | 0:30 | 15:00 |

---

## Slide 1. Title - 0:30

Kính chào thầy và các bạn.

Nhóm em xin trình bày đề tài **Effects of Superspreaders in Spread of Epidemic**, nghĩa là ảnh hưởng của những cá thể siêu lây nhiễm trong quá trình lan truyền dịch bệnh.

Bài tham khảo chính là bài của Fujie và Odagaki, đăng trên *Physica A* năm 2007. Nhóm cài đặt lại mô hình SIR không gian, chạy lại các thí nghiệm Monte Carlo chính, và phân tích vì sao superspreaders có thể làm dịch lan nhanh hơn.

---

## Slide 2. Roadmap - 0:40

Bài trình bày gồm năm phần.

Một là bài toán và lý do cần xét superspreaders. Hai là mô hình SIR không gian với hai cơ chế superspreader. Ba là thuật toán Monte Carlo. Bốn là các kết quả đã chạy lại: percolation, tốc độ lan truyền, đường cong dịch và phân phối ca lây thứ cấp.

Cuối cùng, nhóm thảo luận độ nhạy, tính robust, hạn chế, và kết luận chính từ bài tham khảo.

---

## Slide 3. Problem - 1:00

Trong mô hình SIR cổ điển, dân số được chia thành ba trạng thái: susceptible, infected và recovered. Mô hình trộn đều thường giả sử các cá thể nhiễm bệnh có khả năng lây tương đối giống nhau.

Trong thực tế, nhiều dịch bệnh có tính không đồng nhất mạnh: đa số ca bệnh chỉ lây cho ít người, nhưng một số ít ca tạo ra rất nhiều ca lây thứ cấp. Những cá thể như vậy gọi là **superspreaders**, hay cá thể siêu lây nhiễm.

Bài báo đặt ra một câu hỏi trung tâm: superspreaders quan trọng vì họ lây mạnh hơn ở phạm vi gần, hay vì họ có nhiều kết nối xã hội hơn và có thể truyền bệnh qua các đường tiếp xúc xa hơn?

---

## Slide 4. Modeling Approach - 1:10

Bài báo dùng mô hình SIR trong không gian hai chiều liên tục.

Có \(N\) cá thể đặt cố định trong một hình vuông cạnh \(L\). Trong mô phỏng, \(L=10r_0\), với \(r_0\) là bán kính lây chuẩn của cá thể bình thường.

Mỗi cá thể có một trong ba trạng thái \(S\), \(I\), hoặc \(R\). Mật độ dân số là
\[
\rho=\frac{N}{L^2}.
\]

Điểm khác với SIR trộn đều là xác suất lây phụ thuộc vào khoảng cách. Một cá thể infected lây cho một cá thể susceptible với xác suất \(w(r)\), trong đó \(r\) là khoảng cách tuần hoàn giữa hai cá thể.

Một tỉ lệ \(\lambda\) của dân số được gán là superspreaders, nên \(\lambda\) đo mức độ phổ biến của superspreaders trong quần thể.

---

## Slide 5. Two Superspreader Models - 1:20

Bài báo đề xuất hai cách mô hình hóa superspreaders.

Mô hình thứ nhất là **strong infectiousness**. Superspreader vẫn chỉ lây trong bán kính \(r_0\), nhưng xác suất lây cao hơn. Nếu khoảng cách nhỏ hơn hoặc bằng \(r_0\), xác suất là \(w_0\); ngoài \(r_0\), xác suất bằng 0.

Cách hiểu trực giác là: cá thể này có thể phát tán mầm bệnh mạnh hơn, nhưng ảnh hưởng vẫn mang tính địa phương.

Mô hình thứ hai là **hub model**. Superspreader không nhất thiết lây mạnh hơn ở khoảng cách gần, mà có vùng tiếp xúc hiệu dụng rộng hơn, với bán kính
\[
r_n=\sqrt{6}r_0.
\]

Mô hình này đại diện cho người có nhiều kết nối xã hội, giúp bệnh truyền xa hơn trong không gian.

Hai mô hình được chuẩn hóa để sức lây trung bình có thể so sánh được. Vì vậy, ta đang so sánh cường độ lây cục bộ với độ kết nối xa.

---

## Slide 6. Normal vs Superspreader Infection Kernels - 1:10

Slide này so sánh xác suất lây theo khoảng cách.

Với cá thể nhiễm bệnh bình thường, xác suất lây giảm dần theo khoảng cách:
\[
w_n(r)=w_0(1-r/r_0)^2
\]
khi \(0\le r\le r_0\), và bằng 0 khi ra ngoài \(r_0\).

Trong strong infectiousness, superspreader có xác suất bằng hằng số \(w_0\) trong cùng bán kính \(r_0\). Vì vậy đường màu coral cao hơn trong vùng gần, nhưng không mở rộng ra xa.

Trong hub model, xác suất vẫn giảm theo khoảng cách, nhưng bán kính kéo dài tới \(\sqrt{6}r_0\). Điều này tạo ra các đường lây xa, giúp dịch không chỉ lan qua từng cụm địa phương mà có thể nhảy sang vùng khác.

---

## Slide 7. Basic Reproductive Number and Critical Density - 1:10

Để giải thích ngưỡng lan truyền, bài báo đưa ra đại lượng \(R_0(\lambda)\).

Ở đây, \(R_0(\lambda)\) được tính từ tích phân của các hàm lây theo khoảng cách, gồm nhóm superspreader với tỉ lệ \(\lambda\), và nhóm bình thường với tỉ lệ \(1-\lambda\).

Ý nghĩa trực giác là: \(R_0(\lambda)\) đo số ca lây mới kỳ vọng do một cá thể nhiễm tạo ra trong một bước thời gian, có xét đến mật độ dân số và loại cá thể.

Dịch có thể percolate khi
\[
R_0=R_c.
\]
Từ đó bài báo suy ra mật độ tới hạn:
\[
\rho_c\pi r_0^2=
\frac{R_c}{\lambda+\frac{1-\lambda}{6}}.
\]

Khi \(\lambda\) tăng, mẫu số tăng, nên mật độ tới hạn giảm. Nói cách khác, có thêm superspreaders thì dịch có thể lan qua toàn hệ thống ở mật độ thấp hơn.

---

## Slide 8. Monte Carlo Solving Algorithm - 1:00

Mô hình được giải bằng mô phỏng Monte Carlo.

Ban đầu, nhóm đặt ngẫu nhiên \(N\) cá thể trong miền không gian. Một cá thể là ca nhiễm ban đầu, các cá thể còn lại là susceptible.

Sau đó, một tỉ lệ \(\lambda\) cá thể được gán nhãn superspreader.

Ở mỗi bước thời gian, chỉ những cá thể infected từ đầu bước mới được thử lây. Với từng cặp infected và susceptible, nhóm tính khoảng cách tuần hoàn \(r\), rồi tính xác suất \(w(r)\).

Nếu một số ngẫu nhiên nhỏ hơn \(w(r)\), cá thể susceptible trở thành infected mới. Các ca mới này chỉ bắt đầu lây từ bước sau.

Cuối bước, các cá thể infected hiện tại hồi phục với xác suất \(\gamma\). Trong thí nghiệm tham chiếu, \(\gamma=1\), nên mỗi cá thể chỉ lây trong một bước.

---

## Slide 9. Result 1: Percolation - 1:15

Đây là kết quả đầu tiên: xác suất percolation.

Trục ngang là mật độ chuẩn hóa \(\rho\pi r_0^2\), còn trục dọc là xác suất dịch lan tới phía xa của miền không gian. Với mỗi bộ tham số, nhóm chạy 1000 lần Monte Carlo, theo thiết lập của bài báo.

Khi mật độ thấp, xác suất percolation gần 0. Các chuỗi lây thường tự tắt trước khi lan xa.

Khi mật độ cao, xác suất gần 1. Nghĩa là trong hầu hết mô phỏng, dịch lan qua hệ thống.

Ta thấy khi \(\lambda\) tăng, các đường cong dịch sang trái. Điều này cho thấy mật độ tới hạn giảm khi tỉ lệ superspreader tăng.

Đồ thị critical density xác nhận xu hướng này. Hub model có mật độ tới hạn thấp hơn strong infectiousness model, nghĩa là kết nối xa giúp dịch percolate dễ hơn so với chỉ tăng xác suất lây gần.

---

## Slide 10. Result 2: Propagation Speed - 1:00

Kết quả thứ hai là tốc độ lan truyền.

Đồ thị bên trái biểu diễn front distance \(r_f(t)\), tức khoảng cách xa nhất từ ca nhiễm ban đầu đến các cá thể đã từng nhiễm.

Khi \(\lambda\) tăng, front di chuyển nhanh hơn, nghĩa là superspreaders làm dịch mở rộng nhanh hơn.

Đồ thị bên phải so sánh vận tốc ước lượng giữa hai mô hình. Với \(\lambda>0\), hub model có tốc độ lớn hơn.

Lý do là hub superspreaders tạo ra đường lây tầm xa. Thay vì chỉ lan dần qua các cụm gần nhau, dịch có thể nhảy sang khu vực xa và tạo vùng nhiễm mới sớm hơn.

---

## Slide 11. Result 3: Epidemic Curves - 0:55

Kết quả thứ ba là đường cong dịch, tức số ca nhiễm mới theo từng bước thời gian.

Khi không có superspreaders, đường cong rộng hơn và đỉnh xuất hiện muộn hơn, tương ứng với lan truyền chậm hơn.

Khi có superspreaders, dịch tăng nhanh hơn và đạt đỉnh sớm hơn.

Trong các đường cong này, hub model có pha tăng nhanh nhất. Điều này phù hợp với kết quả tốc độ: kết nối xa không chỉ làm dịch lan xa hơn, mà còn làm bùng phát sớm hơn.

---

## Slide 12. Result 4: Secondary Infections - 0:55

Kết quả tiếp theo là phân phối số ca lây thứ cấp.

Trong mạng lây nhiễm, mỗi cạnh có hướng biểu diễn việc một cá thể lây cho cá thể khác. Số cạnh đi ra chính là số ca lây thứ cấp của cá thể đó.

Khi không có superspreaders, phần lớn cá thể gây ra rất ít ca lây, và phân phối có đuôi ngắn.

Khi có superspreaders, đa số cá thể vẫn không lây nhiều, nhưng xuất hiện một nhóm nhỏ tạo ra rất nhiều ca lây thứ cấp. Vì vậy phân phối có đuôi dài.

Đuôi dài này là dấu hiệu thống kê quan trọng của hiện tượng superspreading.

---

## Slide 13. Comparison with SARS Data - 0:55

Bài báo cũng so sánh mô hình với dữ liệu SARS tại Singapore năm 2003.

Trong dữ liệu thực tế, có một số ca gây ra rất nhiều ca lây thứ cấp, ví dụ 12, 21, 23 và 40 ca. Đây là kiểu dữ liệu đuôi dài, giống mô phỏng khi có superspreaders.

Bài báo nhận thấy đường cong dịch phù hợp tốt nhất về mặt định tính khi dùng hub model với
\[
N=477,\quad \rho\pi r_0^2=15,\quad \lambda=0.4,\quad \gamma=1.
\]

Vì vậy, bài tham khảo kết luận rằng superspreading trong SARS được giải thích tốt hơn bởi nhiều kết nối xã hội, thay vì chỉ bởi khả năng lây cục bộ mạnh hơn.

---

## Slide 14. Sensitivity and Robustness - 0:55

Mô hình nhạy với một số tham số chính.

Thứ nhất là tỉ lệ superspreader \(\lambda\). Khi \(\lambda\) tăng, mật độ tới hạn giảm và tốc độ lan truyền tăng.

Thứ hai là cơ chế superspreading. Hai mô hình có sức lây trung bình được chuẩn hóa, nhưng hub model vẫn lan nhanh hơn vì tạo ra liên kết xa.

Thứ ba là mật độ dân số. Ở mật độ thấp, chuỗi lây dễ bị dập tắt. Ở mật độ cao, mạng tiếp xúc đủ dày để dịch percolate.

Cuối cùng là xác suất hồi phục \(\gamma\). Nếu \(\gamma\) thấp hơn, thời gian infectious dài hơn, nên có nhiều cơ hội lây.

Kết luận robust là: superspreaders làm giảm ngưỡng lan truyền, tăng tốc độ lan, làm đỉnh dịch đến sớm hơn, và tạo phân phối ca lây thứ cấp có đuôi dài.

---

## Slide 15. Limitations - 0:35

Mô hình cũng có một số hạn chế.

Thứ nhất, các cá thể được đặt cố định, trong khi thực tế con người di chuyển và mạng tiếp xúc thay đổi theo thời gian.

Thứ hai, nhãn superspreader được gán ngẫu nhiên, nhưng thực tế nó có thể phụ thuộc vào nghề nghiệp, hành vi, môi trường, sinh học, hoặc sự kiện đông người.

Vì vậy, mô hình không phải công cụ dự báo đầy đủ cho mọi dịch bệnh. Giá trị chính là so sánh hai cơ chế superspreading trong một khung toán học có kiểm soát.

---

## Slide 16. Conclusion - 0:30

Tóm lại, superspreaders làm thay đổi đáng kể động lực lan truyền dịch.

Khi tỉ lệ superspreader tăng, mật độ tới hạn giảm, tốc độ lan truyền tăng, đỉnh dịch đến sớm hơn, và phân phối số ca lây thứ cấp có đuôi dài.

Giữa hai cơ chế, hub model lan nhanh hơn và phù hợp hơn với so sánh SARS trong bài báo.

Thông điệp chính là: chỉ dùng hệ số lây trung bình là chưa đủ. Mô hình dịch bệnh cần xét đến sự không đồng nhất trong kết nối giữa các cá thể.

Em xin cảm ơn thầy và các bạn đã lắng nghe.
