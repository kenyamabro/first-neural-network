# first-neural-network
I created a neural network without using specific libraries dedicated to neural networks, implementing every neural network calculus directly in code, including initialization, activation, backpropagation and minimization by negative gradient's descent.
The hidden layers number and the neurons number inside each of them can be easily modified by changing a QList of int : neuronsNum (mainwindow.h l.25). The default list is {728, 16, 16, 10}.
This was not created watching a programming tutorial telling me exactly how to do it. I studied and learned the algorithm first, and then translated it to code.

One example of cost functions and accuracy of the neural network over 600 mini-batches of 100 data points:  
0 # Cost: 0.552054 Accuracy: 2 / 100  
1 # Cost: 0.512597 Accuracy: 3 / 100  
2 # Cost: 0.515335 Accuracy: 5 / 100  
3 # Cost: 0.471051 Accuracy: 8 / 100  
4 # Cost: 0.461714 Accuracy: 4 / 100  
5 # Cost: 0.425692 Accuracy: 2 / 100  
6 # Cost: 0.417987 Accuracy: 8 / 100  
7 # Cost: 0.402937 Accuracy: 8 / 100  
8 # Cost: 0.385421 Accuracy: 10 / 100  
9 # Cost: 0.369186 Accuracy: 15 / 100  
10 # Cost: 0.363857 Accuracy: 13 / 100  
11 # Cost: 0.348435 Accuracy: 14 / 100  
12 # Cost: 0.318671 Accuracy: 18 / 100  
13 # Cost: 0.301252 Accuracy: 26 / 100  
14 # Cost: 0.288691 Accuracy: 27 / 100  
15 # Cost: 0.273903 Accuracy: 27 / 100  
16 # Cost: 0.269628 Accuracy: 37 / 100  
17 # Cost: 0.247562 Accuracy: 36 / 100  
18 # Cost: 0.239116 Accuracy: 44 / 100  
19 # Cost: 0.227613 Accuracy: 46 / 100  
20 # Cost: 0.215799 Accuracy: 58 / 100  
21 # Cost: 0.2002 Accuracy: 60 / 100  
22 # Cost: 0.186519 Accuracy: 65 / 100  
23 # Cost: 0.180578 Accuracy: 68 / 100  
24 # Cost: 0.175339 Accuracy: 74 / 100  
25 # Cost: 0.162897 Accuracy: 78 / 100  
26 # Cost: 0.16839 Accuracy: 74 / 100  
27 # Cost: 0.162063 Accuracy: 77 / 100  
28 # Cost: 0.158718 Accuracy: 85 / 100  
29 # Cost: 0.150408 Accuracy: 79 / 100  
30 # Cost: 0.145416 Accuracy: 87 / 100  
31 # Cost: 0.139761 Accuracy: 87 / 100  
32 # Cost: 0.14375 Accuracy: 91 / 100  
33 # Cost: 0.142965 Accuracy: 92 / 100  
34 # Cost: 0.125596 Accuracy: 84 / 100  
35 # Cost: 0.133883 Accuracy: 89 / 100  
36 # Cost: 0.130528 Accuracy: 91 / 100  
37 # Cost: 0.131391 Accuracy: 92 / 100  
38 # Cost: 0.125724 Accuracy: 93 / 100  
39 # Cost: 0.116306 Accuracy: 92 / 100  
40 # Cost: 0.120587 Accuracy: 94 / 100  
41 # Cost: 0.127408 Accuracy: 96 / 100  
42 # Cost: 0.11964 Accuracy: 94 / 100  
43 # Cost: 0.12167 Accuracy: 97 / 100  
44 # Cost: 0.118038 Accuracy: 96 / 100  
45 # Cost: 0.114365 Accuracy: 98 / 100  
46 # Cost: 0.116208 Accuracy: 97 / 100  
47 # Cost: 0.123607 Accuracy: 96 / 100  
48 # Cost: 0.117676 Accuracy: 97 / 100  
49 # Cost: 0.115246 Accuracy: 99 / 100  
50 # Cost: 0.114062 Accuracy: 100 / 100  
51 # Cost: 0.115291 Accuracy: 99 / 100  
52 # Cost: 0.113853 Accuracy: 99 / 100  
53 # Cost: 0.116181 Accuracy: 100 / 100  
54 # Cost: 0.103955 Accuracy: 100 / 100  
55 # Cost: 0.11577 Accuracy: 97 / 100  
56 # Cost: 0.113288 Accuracy: 95 / 100  
57 # Cost: 0.106712 Accuracy: 98 / 100  
58 # Cost: 0.114261 Accuracy: 100 / 100  
59 # Cost: 0.106711 Accuracy: 97 / 100  
60 # Cost: 0.104872 Accuracy: 98 / 100  
61 # Cost: 0.10381 Accuracy: 100 / 100  
62 # Cost: 0.102861 Accuracy: 100 / 100  
63 # Cost: 0.105434 Accuracy: 98 / 100  
64 # Cost: 0.101175 Accuracy: 100 / 100  
65 # Cost: 0.107029 Accuracy: 99 / 100  
66 # Cost: 0.103681 Accuracy: 99 / 100  
67 # Cost: 0.105533 Accuracy: 100 / 100  
68 # Cost: 0.103212 Accuracy: 98 / 100  
69 # Cost: 0.106075 Accuracy: 99 / 100  
70 # Cost: 0.105124 Accuracy: 100 / 100  
71 # Cost: 0.108272 Accuracy: 100 / 100  
72 # Cost: 0.104722 Accuracy: 98 / 100  
73 # Cost: 0.103481 Accuracy: 99 / 100  
74 # Cost: 0.0993988 Accuracy: 98 / 100  
75 # Cost: 0.100318 Accuracy: 100 / 100  
76 # Cost: 0.101438 Accuracy: 99 / 100  
77 # Cost: 0.101492 Accuracy: 99 / 100  
78 # Cost: 0.105168 Accuracy: 99 / 100  
79 # Cost: 0.100072 Accuracy: 99 / 100  
80 # Cost: 0.102909 Accuracy: 100 / 100  
81 # Cost: 0.103731 Accuracy: 100 / 100  
82 # Cost: 0.104922 Accuracy: 100 / 100  
83 # Cost: 0.0962725 Accuracy: 100 / 100  
84 # Cost: 0.102385 Accuracy: 99 / 100  
85 # Cost: 0.0986282 Accuracy: 100 / 100  
86 # Cost: 0.100113 Accuracy: 97 / 100  
87 # Cost: 0.0992205 Accuracy: 100 / 100  
88 # Cost: 0.10289 Accuracy: 99 / 100  
89 # Cost: 0.0981689 Accuracy: 100 / 100  
90 # Cost: 0.0977622 Accuracy: 99 / 100  
91 # Cost: 0.100974 Accuracy: 100 / 100  
92 # Cost: 0.101101 Accuracy: 98 / 100  
93 # Cost: 0.098439 Accuracy: 100 / 100  
94 # Cost: 0.101685 Accuracy: 100 / 100  
95 # Cost: 0.101164 Accuracy: 100 / 100  
96 # Cost: 0.102344 Accuracy: 100 / 100  
97 # Cost: 0.0975371 Accuracy: 99 / 100  
98 # Cost: 0.102103 Accuracy: 99 / 100  
99 # Cost: 0.0975502 Accuracy: 100 / 100  
100 # Cost: 0.0984319 Accuracy: 100 / 100  
101 # Cost: 0.0991831 Accuracy: 100 / 100  
102 # Cost: 0.0952514 Accuracy: 100 / 100  
103 # Cost: 0.0996027 Accuracy: 100 / 100  
104 # Cost: 0.0966838 Accuracy: 99 / 100  
105 # Cost: 0.0975584 Accuracy: 100 / 100  
106 # Cost: 0.100925 Accuracy: 100 / 100  
107 # Cost: 0.101023 Accuracy: 99 / 100  
108 # Cost: 0.0973233 Accuracy: 100 / 100  
109 # Cost: 0.09944 Accuracy: 100 / 100  
110 # Cost: 0.10079 Accuracy: 100 / 100  
111 # Cost: 0.0959239 Accuracy: 100 / 100  
112 # Cost: 0.0974513 Accuracy: 100 / 100  
113 # Cost: 0.0960868 Accuracy: 100 / 100  
114 # Cost: 0.0988623 Accuracy: 100 / 100  
115 # Cost: 0.0955449 Accuracy: 100 / 100  
116 # Cost: 0.0982298 Accuracy: 99 / 100  
117 # Cost: 0.0965341 Accuracy: 100 / 100  
118 # Cost: 0.094539 Accuracy: 100 / 100  
119 # Cost: 0.0948142 Accuracy: 100 / 100  
120 # Cost: 0.0985983 Accuracy: 100 / 100  
121 # Cost: 0.0987321 Accuracy: 100 / 100  
122 # Cost: 0.0968067 Accuracy: 100 / 100  
123 # Cost: 0.0961533 Accuracy: 100 / 100  
124 # Cost: 0.0970441 Accuracy: 100 / 100  
125 # Cost: 0.098933 Accuracy: 99 / 100  
126 # Cost: 0.0945928 Accuracy: 100 / 100  
127 # Cost: 0.0945421 Accuracy: 100 / 100  
128 # Cost: 0.0956383 Accuracy: 99 / 100  
129 # Cost: 0.0961724 Accuracy: 100 / 100  
130 # Cost: 0.0938864 Accuracy: 100 / 100  
131 # Cost: 0.0982234 Accuracy: 100 / 100  
132 # Cost: 0.0950278 Accuracy: 99 / 100  
133 # Cost: 0.0965783 Accuracy: 100 / 100  
134 # Cost: 0.0958013 Accuracy: 100 / 100  
135 # Cost: 0.0947822 Accuracy: 100 / 100  
136 # Cost: 0.0974646 Accuracy: 100 / 100  
137 # Cost: 0.0939298 Accuracy: 100 / 100  
138 # Cost: 0.0992718 Accuracy: 100 / 100  
139 # Cost: 0.094166 Accuracy: 99 / 100  
140 # Cost: 0.0982923 Accuracy: 100 / 100  
141 # Cost: 0.0995132 Accuracy: 100 / 100  
142 # Cost: 0.097286 Accuracy: 100 / 100  
143 # Cost: 0.0936934 Accuracy: 100 / 100  
144 # Cost: 0.0960766 Accuracy: 100 / 100  
145 # Cost: 0.0949271 Accuracy: 100 / 100  
146 # Cost: 0.10102 Accuracy: 100 / 100  
147 # Cost: 0.100776 Accuracy: 100 / 100  
148 # Cost: 0.0986237 Accuracy: 100 / 100  
149 # Cost: 0.0967117 Accuracy: 100 / 100  
150 # Cost: 0.0958087 Accuracy: 100 / 100  
151 # Cost: 0.0962483 Accuracy: 100 / 100  
152 # Cost: 0.0951481 Accuracy: 100 / 100  
153 # Cost: 0.0950814 Accuracy: 100 / 100  
154 # Cost: 0.0949059 Accuracy: 100 / 100  
155 # Cost: 0.0966791 Accuracy: 100 / 100  
156 # Cost: 0.0932994 Accuracy: 100 / 100  
157 # Cost: 0.0940334 Accuracy: 99 / 100  
158 # Cost: 0.0948287 Accuracy: 100 / 100  
159 # Cost: 0.0947185 Accuracy: 100 / 100  
160 # Cost: 0.0932727 Accuracy: 99 / 100  
161 # Cost: 0.0943154 Accuracy: 100 / 100  
162 # Cost: 0.0956011 Accuracy: 100 / 100  
163 # Cost: 0.093998 Accuracy: 100 / 100  
164 # Cost: 0.0949153 Accuracy: 100 / 100  
165 # Cost: 0.0955338 Accuracy: 100 / 100  
166 # Cost: 0.0947971 Accuracy: 100 / 100  
167 # Cost: 0.0940615 Accuracy: 100 / 100  
168 # Cost: 0.0950382 Accuracy: 100 / 100  
169 # Cost: 0.0932963 Accuracy: 100 / 100  
170 # Cost: 0.0956159 Accuracy: 100 / 100  
171 # Cost: 0.0938617 Accuracy: 100 / 100  
172 # Cost: 0.0953673 Accuracy: 100 / 100  
173 # Cost: 0.0938662 Accuracy: 100 / 100  
174 # Cost: 0.0945579 Accuracy: 100 / 100  
175 # Cost: 0.0927081 Accuracy: 99 / 100  
176 # Cost: 0.0929839 Accuracy: 100 / 100  
177 # Cost: 0.0941245 Accuracy: 100 / 100  
178 # Cost: 0.0949218 Accuracy: 100 / 100  
179 # Cost: 0.0943223 Accuracy: 100 / 100  
180 # Cost: 0.0931609 Accuracy: 100 / 100  
181 # Cost: 0.0927063 Accuracy: 100 / 100  
182 # Cost: 0.0933603 Accuracy: 100 / 100  
183 # Cost: 0.0935559 Accuracy: 99 / 100  
184 # Cost: 0.0926882 Accuracy: 100 / 100  
185 # Cost: 0.0930985 Accuracy: 100 / 100  
186 # Cost: 0.096063 Accuracy: 100 / 100  
187 # Cost: 0.096288 Accuracy: 100 / 100  
188 # Cost: 0.0960121 Accuracy: 100 / 100  
189 # Cost: 0.0959804 Accuracy: 100 / 100  
190 # Cost: 0.0911625 Accuracy: 100 / 100  
191 # Cost: 0.0933697 Accuracy: 100 / 100  
192 # Cost: 0.0931726 Accuracy: 100 / 100  
193 # Cost: 0.0939357 Accuracy: 100 / 100  
194 # Cost: 0.093211 Accuracy: 100 / 100  
195 # Cost: 0.0929183 Accuracy: 100 / 100  
196 # Cost: 0.0910541 Accuracy: 100 / 100  
197 # Cost: 0.0910195 Accuracy: 100 / 100  
198 # Cost: 0.0924601 Accuracy: 99 / 100  
199 # Cost: 0.0945409 Accuracy: 100 / 100  
200 # Cost: 0.0939986 Accuracy: 100 / 100  
201 # Cost: 0.0955166 Accuracy: 100 / 100  
202 # Cost: 0.0952012 Accuracy: 100 / 100  
203 # Cost: 0.0947571 Accuracy: 100 / 100  
204 # Cost: 0.0922657 Accuracy: 99 / 100  
205 # Cost: 0.0919489 Accuracy: 100 / 100  
206 # Cost: 0.0947238 Accuracy: 100 / 100  
207 # Cost: 0.0945044 Accuracy: 100 / 100  
208 # Cost: 0.095441 Accuracy: 100 / 100  
209 # Cost: 0.0935887 Accuracy: 100 / 100  
210 # Cost: 0.0952837 Accuracy: 99 / 100  
211 # Cost: 0.0936628 Accuracy: 100 / 100  
212 # Cost: 0.0929347 Accuracy: 100 / 100  
213 # Cost: 0.0932843 Accuracy: 100 / 100  
214 # Cost: 0.092467 Accuracy: 100 / 100  
215 # Cost: 0.0913611 Accuracy: 100 / 100  
216 # Cost: 0.0925697 Accuracy: 100 / 100  
217 # Cost: 0.0936631 Accuracy: 100 / 100  
218 # Cost: 0.0909279 Accuracy: 100 / 100  
219 # Cost: 0.0933499 Accuracy: 100 / 100  
220 # Cost: 0.0957535 Accuracy: 100 / 100  
221 # Cost: 0.099624 Accuracy: 100 / 100  
222 # Cost: 0.0950297 Accuracy: 100 / 100  
223 # Cost: 0.0943868 Accuracy: 100 / 100  
224 # Cost: 0.0934574 Accuracy: 100 / 100  
225 # Cost: 0.0921744 Accuracy: 100 / 100  
226 # Cost: 0.0927374 Accuracy: 100 / 100  
227 # Cost: 0.0937923 Accuracy: 100 / 100  
228 # Cost: 0.0943231 Accuracy: 100 / 100  
229 # Cost: 0.0940823 Accuracy: 100 / 100  
230 # Cost: 0.0930722 Accuracy: 100 / 100  
231 # Cost: 0.0927797 Accuracy: 100 / 100  
232 # Cost: 0.0915568 Accuracy: 100 / 100  
233 # Cost: 0.0946652 Accuracy: 100 / 100  
234 # Cost: 0.0931052 Accuracy: 100 / 100  
235 # Cost: 0.0922153 Accuracy: 100 / 100  
236 # Cost: 0.0938279 Accuracy: 100 / 100  
237 # Cost: 0.0930289 Accuracy: 100 / 100  
238 # Cost: 0.0941024 Accuracy: 100 / 100  
239 # Cost: 0.0924642 Accuracy: 100 / 100  
240 # Cost: 0.0934028 Accuracy: 100 / 100  
241 # Cost: 0.0931362 Accuracy: 100 / 100  
242 # Cost: 0.0927677 Accuracy: 100 / 100  
243 # Cost: 0.0937585 Accuracy: 100 / 100  
244 # Cost: 0.0950323 Accuracy: 100 / 100  
245 # Cost: 0.0932097 Accuracy: 100 / 100  
246 # Cost: 0.093495 Accuracy: 100 / 100  
247 # Cost: 0.0923173 Accuracy: 100 / 100  
248 # Cost: 0.0936977 Accuracy: 100 / 100  
249 # Cost: 0.094787 Accuracy: 100 / 100  
250 # Cost: 0.0909861 Accuracy: 100 / 100  
251 # Cost: 0.0921062 Accuracy: 100 / 100  
252 # Cost: 0.0925689 Accuracy: 100 / 100  
253 # Cost: 0.0932942 Accuracy: 100 / 100  
254 # Cost: 0.0909792 Accuracy: 100 / 100  
255 # Cost: 0.0925184 Accuracy: 100 / 100  
256 # Cost: 0.0929407 Accuracy: 100 / 100  
257 # Cost: 0.0937473 Accuracy: 100 / 100  
258 # Cost: 0.0927667 Accuracy: 100 / 100  
259 # Cost: 0.0924958 Accuracy: 100 / 100  
260 # Cost: 0.0934356 Accuracy: 100 / 100  
261 # Cost: 0.0933577 Accuracy: 100 / 100  
262 # Cost: 0.094521 Accuracy: 100 / 100  
263 # Cost: 0.0929299 Accuracy: 100 / 100  
264 # Cost: 0.0923152 Accuracy: 100 / 100  
265 # Cost: 0.0921639 Accuracy: 99 / 100  
266 # Cost: 0.0912231 Accuracy: 100 / 100  
267 # Cost: 0.092077 Accuracy: 100 / 100  
268 # Cost: 0.0923068 Accuracy: 100 / 100  
269 # Cost: 0.0929546 Accuracy: 100 / 100  
270 # Cost: 0.0915527 Accuracy: 100 / 100  
271 # Cost: 0.0914652 Accuracy: 100 / 100  
272 # Cost: 0.0927922 Accuracy: 100 / 100  
273 # Cost: 0.0936841 Accuracy: 100 / 100  
274 # Cost: 0.0928791 Accuracy: 100 / 100  
275 # Cost: 0.0944758 Accuracy: 100 / 100  
276 # Cost: 0.0939336 Accuracy: 100 / 100  
277 # Cost: 0.0932915 Accuracy: 100 / 100  
278 # Cost: 0.0925778 Accuracy: 100 / 100  
279 # Cost: 0.0922812 Accuracy: 100 / 100  
280 # Cost: 0.0922799 Accuracy: 100 / 100  
281 # Cost: 0.0936997 Accuracy: 100 / 100  
282 # Cost: 0.092481 Accuracy: 100 / 100  
283 # Cost: 0.0921867 Accuracy: 100 / 100  
284 # Cost: 0.0922285 Accuracy: 99 / 100  
285 # Cost: 0.090158 Accuracy: 100 / 100  
286 # Cost: 0.0932035 Accuracy: 100 / 100  
287 # Cost: 0.0923246 Accuracy: 100 / 100  
288 # Cost: 0.093564 Accuracy: 100 / 100  
289 # Cost: 0.0905358 Accuracy: 100 / 100  
290 # Cost: 0.0924084 Accuracy: 100 / 100  
291 # Cost: 0.0922593 Accuracy: 100 / 100  
292 # Cost: 0.0941528 Accuracy: 100 / 100  
293 # Cost: 0.0924785 Accuracy: 100 / 100  
294 # Cost: 0.0914942 Accuracy: 100 / 100  
295 # Cost: 0.0929973 Accuracy: 100 / 100  
296 # Cost: 0.0923415 Accuracy: 100 / 100  
297 # Cost: 0.0923641 Accuracy: 99 / 100  
298 # Cost: 0.0929592 Accuracy: 100 / 100  
299 # Cost: 0.0909828 Accuracy: 100 / 100  
300 # Cost: 0.0926121 Accuracy: 100 / 100  
301 # Cost: 0.0943377 Accuracy: 100 / 100  
302 # Cost: 0.0940894 Accuracy: 100 / 100  
303 # Cost: 0.0925517 Accuracy: 100 / 100  
304 # Cost: 0.0917239 Accuracy: 100 / 100  
305 # Cost: 0.0921885 Accuracy: 100 / 100  
306 # Cost: 0.0939862 Accuracy: 100 / 100  
307 # Cost: 0.0918121 Accuracy: 100 / 100  
308 # Cost: 0.092646 Accuracy: 100 / 100  
309 # Cost: 0.0924023 Accuracy: 100 / 100  
310 # Cost: 0.0927951 Accuracy: 100 / 100  
311 # Cost: 0.0941242 Accuracy: 100 / 100  
312 # Cost: 0.0932781 Accuracy: 100 / 100  
313 # Cost: 0.0917964 Accuracy: 100 / 100  
314 # Cost: 0.0897074 Accuracy: 100 / 100  
315 # Cost: 0.0916499 Accuracy: 100 / 100  
316 # Cost: 0.0934194 Accuracy: 100 / 100  
317 # Cost: 0.0928403 Accuracy: 100 / 100  
318 # Cost: 0.0918651 Accuracy: 100 / 100  
319 # Cost: 0.0911541 Accuracy: 100 / 100  
320 # Cost: 0.0915117 Accuracy: 100 / 100  
321 # Cost: 0.0930187 Accuracy: 100 / 100  
322 # Cost: 0.0914602 Accuracy: 100 / 100  
323 # Cost: 0.0929305 Accuracy: 100 / 100  
324 # Cost: 0.0921764 Accuracy: 100 / 100  
325 # Cost: 0.0924056 Accuracy: 100 / 100  
326 # Cost: 0.0892929 Accuracy: 100 / 100  
327 # Cost: 0.0921771 Accuracy: 100 / 100  
328 # Cost: 0.0910168 Accuracy: 100 / 100  
329 # Cost: 0.0926709 Accuracy: 100 / 100  
330 # Cost: 0.0923217 Accuracy: 100 / 100  
331 # Cost: 0.0929543 Accuracy: 100 / 100  
332 # Cost: 0.0922132 Accuracy: 100 / 100  
333 # Cost: 0.0938021 Accuracy: 100 / 100  
334 # Cost: 0.0916801 Accuracy: 100 / 100  
335 # Cost: 0.0918678 Accuracy: 100 / 100  
336 # Cost: 0.0925919 Accuracy: 100 / 100  
337 # Cost: 0.0919242 Accuracy: 100 / 100  
338 # Cost: 0.091806 Accuracy: 100 / 100  
339 # Cost: 0.0908386 Accuracy: 100 / 100  
340 # Cost: 0.0937693 Accuracy: 100 / 100  
341 # Cost: 0.0915588 Accuracy: 100 / 100  
342 # Cost: 0.0930938 Accuracy: 100 / 100  
343 # Cost: 0.0933702 Accuracy: 100 / 100  
344 # Cost: 0.0924887 Accuracy: 100 / 100  
345 # Cost: 0.090833 Accuracy: 100 / 100  
346 # Cost: 0.0915805 Accuracy: 100 / 100  
347 # Cost: 0.0921436 Accuracy: 100 / 100  
348 # Cost: 0.0919935 Accuracy: 100 / 100  
349 # Cost: 0.0905497 Accuracy: 100 / 100  
350 # Cost: 0.0947132 Accuracy: 100 / 100  
351 # Cost: 0.0920357 Accuracy: 100 / 100  
352 # Cost: 0.0941929 Accuracy: 100 / 100  
353 # Cost: 0.0913324 Accuracy: 100 / 100  
354 # Cost: 0.0924272 Accuracy: 100 / 100  
355 # Cost: 0.0922157 Accuracy: 100 / 100  
356 # Cost: 0.0930307 Accuracy: 100 / 100  
357 # Cost: 0.0924387 Accuracy: 99 / 100  
358 # Cost: 0.0913008 Accuracy: 100 / 100  
359 # Cost: 0.0936357 Accuracy: 100 / 100  
360 # Cost: 0.0940954 Accuracy: 100 / 100  
361 # Cost: 0.0915893 Accuracy: 100 / 100  
362 # Cost: 0.092917 Accuracy: 100 / 100  
363 # Cost: 0.0908537 Accuracy: 100 / 100  
364 # Cost: 0.0928043 Accuracy: 100 / 100  
365 # Cost: 0.0911311 Accuracy: 100 / 100  
366 # Cost: 0.0928072 Accuracy: 100 / 100  
367 # Cost: 0.0918228 Accuracy: 100 / 100  
368 # Cost: 0.0919424 Accuracy: 100 / 100  
369 # Cost: 0.0923265 Accuracy: 100 / 100  
370 # Cost: 0.0916807 Accuracy: 100 / 100  
371 # Cost: 0.0930155 Accuracy: 100 / 100  
372 # Cost: 0.091603 Accuracy: 100 / 100  
373 # Cost: 0.0919635 Accuracy: 100 / 100  
374 # Cost: 0.0920733 Accuracy: 100 / 100  
375 # Cost: 0.0921222 Accuracy: 100 / 100  
376 # Cost: 0.0935605 Accuracy: 100 / 100  
377 # Cost: 0.0918003 Accuracy: 99 / 100  
378 # Cost: 0.0933082 Accuracy: 100 / 100  
379 # Cost: 0.0917168 Accuracy: 100 / 100  
380 # Cost: 0.0917295 Accuracy: 100 / 100  
381 # Cost: 0.0932446 Accuracy: 100 / 100  
382 # Cost: 0.0925415 Accuracy: 100 / 100  
383 # Cost: 0.0907085 Accuracy: 100 / 100  
384 # Cost: 0.092711 Accuracy: 100 / 100  
385 # Cost: 0.0910615 Accuracy: 100 / 100  
386 # Cost: 0.0915882 Accuracy: 100 / 100  
387 # Cost: 0.0913599 Accuracy: 100 / 100  
388 # Cost: 0.0940317 Accuracy: 100 / 100  
389 # Cost: 0.0918487 Accuracy: 100 / 100  
390 # Cost: 0.0919495 Accuracy: 100 / 100  
391 # Cost: 0.0921599 Accuracy: 100 / 100  
392 # Cost: 0.0917205 Accuracy: 100 / 100  
393 # Cost: 0.0922438 Accuracy: 100 / 100  
394 # Cost: 0.0935583 Accuracy: 100 / 100  
395 # Cost: 0.0908023 Accuracy: 100 / 100  
396 # Cost: 0.0924139 Accuracy: 100 / 100  
397 # Cost: 0.0920799 Accuracy: 100 / 100  
398 # Cost: 0.0923584 Accuracy: 100 / 100  
399 # Cost: 0.0923966 Accuracy: 100 / 100  
400 # Cost: 0.0919211 Accuracy: 100 / 100  
401 # Cost: 0.0931346 Accuracy: 100 / 100  
402 # Cost: 0.0904597 Accuracy: 100 / 100  
403 # Cost: 0.090775 Accuracy: 100 / 100  
404 # Cost: 0.0926641 Accuracy: 100 / 100  
405 # Cost: 0.0916834 Accuracy: 100 / 100  
406 # Cost: 0.0926304 Accuracy: 100 / 100  
407 # Cost: 0.092819 Accuracy: 100 / 100  
408 # Cost: 0.0925025 Accuracy: 100 / 100  
409 # Cost: 0.0929094 Accuracy: 100 / 100  
410 # Cost: 0.0899808 Accuracy: 100 / 100  
411 # Cost: 0.0917935 Accuracy: 100 / 100  
412 # Cost: 0.0922478 Accuracy: 100 / 100  
413 # Cost: 0.0927702 Accuracy: 100 / 100  
414 # Cost: 0.0923474 Accuracy: 100 / 100  
415 # Cost: 0.0925393 Accuracy: 100 / 100  
416 # Cost: 0.0907644 Accuracy: 100 / 100  
417 # Cost: 0.0930821 Accuracy: 100 / 100  
418 # Cost: 0.0930708 Accuracy: 100 / 100  
419 # Cost: 0.091691 Accuracy: 100 / 100  
420 # Cost: 0.0911609 Accuracy: 100 / 100  
421 # Cost: 0.0921572 Accuracy: 100 / 100  
422 # Cost: 0.0935692 Accuracy: 100 / 100  
423 # Cost: 0.0923371 Accuracy: 100 / 100  
424 # Cost: 0.0923922 Accuracy: 100 / 100  
425 # Cost: 0.0894537 Accuracy: 100 / 100  
426 # Cost: 0.0932376 Accuracy: 100 / 100  
427 # Cost: 0.0914815 Accuracy: 100 / 100  
428 # Cost: 0.0925503 Accuracy: 100 / 100  
429 # Cost: 0.0921781 Accuracy: 100 / 100  
430 # Cost: 0.092905 Accuracy: 100 / 100  
431 # Cost: 0.091309 Accuracy: 100 / 100  
432 # Cost: 0.0919566 Accuracy: 100 / 100  
433 # Cost: 0.089761 Accuracy: 100 / 100  
434 # Cost: 0.0926162 Accuracy: 100 / 100  
435 # Cost: 0.0914519 Accuracy: 100 / 100  
436 # Cost: 0.0917558 Accuracy: 100 / 100  
437 # Cost: 0.0913798 Accuracy: 100 / 100  
438 # Cost: 0.09218 Accuracy: 100 / 100  
439 # Cost: 0.0896542 Accuracy: 100 / 100  
440 # Cost: 0.0899247 Accuracy: 100 / 100  
441 # Cost: 0.0919272 Accuracy: 100 / 100  
442 # Cost: 0.0937325 Accuracy: 100 / 100  
443 # Cost: 0.0917892 Accuracy: 100 / 100  
444 # Cost: 0.0916935 Accuracy: 100 / 100  
445 # Cost: 0.0911989 Accuracy: 100 / 100  
446 # Cost: 0.0914202 Accuracy: 100 / 100  
447 # Cost: 0.0909387 Accuracy: 100 / 100  
448 # Cost: 0.0915417 Accuracy: 100 / 100  
449 # Cost: 0.091284 Accuracy: 100 / 100  
450 # Cost: 0.0919894 Accuracy: 100 / 100  
451 # Cost: 0.0908823 Accuracy: 100 / 100  
452 # Cost: 0.0920097 Accuracy: 100 / 100  
453 # Cost: 0.0917642 Accuracy: 100 / 100  
454 # Cost: 0.0914283 Accuracy: 100 / 100  
455 # Cost: 0.0910341 Accuracy: 100 / 100  
456 # Cost: 0.0901705 Accuracy: 100 / 100  
457 # Cost: 0.092091 Accuracy: 100 / 100  
458 # Cost: 0.0902878 Accuracy: 100 / 100  
459 # Cost: 0.0930237 Accuracy: 100 / 100  
460 # Cost: 0.0920113 Accuracy: 100 / 100  
461 # Cost: 0.0920397 Accuracy: 100 / 100  
462 # Cost: 0.0934995 Accuracy: 100 / 100  
463 # Cost: 0.0906945 Accuracy: 100 / 100  
464 # Cost: 0.0928668 Accuracy: 100 / 100  
465 # Cost: 0.0925489 Accuracy: 100 / 100  
466 # Cost: 0.091729 Accuracy: 100 / 100  
467 # Cost: 0.0915872 Accuracy: 100 / 100  
468 # Cost: 0.092753 Accuracy: 100 / 100  
469 # Cost: 0.0911482 Accuracy: 100 / 100  
470 # Cost: 0.0920092 Accuracy: 100 / 100  
471 # Cost: 0.0925858 Accuracy: 100 / 100  
472 # Cost: 0.0916142 Accuracy: 100 / 100  
473 # Cost: 0.0912687 Accuracy: 100 / 100  
474 # Cost: 0.0927037 Accuracy: 100 / 100  
475 # Cost: 0.0896401 Accuracy: 100 / 100  
476 # Cost: 0.0935595 Accuracy: 100 / 100  
477 # Cost: 0.091473 Accuracy: 100 / 100  
478 # Cost: 0.091781 Accuracy: 100 / 100  
479 # Cost: 0.0931942 Accuracy: 100 / 100  
480 # Cost: 0.0910248 Accuracy: 100 / 100  
481 # Cost: 0.091352 Accuracy: 100 / 100  
482 # Cost: 0.0916611 Accuracy: 100 / 100  
483 # Cost: 0.0919258 Accuracy: 100 / 100  
484 # Cost: 0.0927621 Accuracy: 100 / 100  
485 # Cost: 0.0913191 Accuracy: 100 / 100  
486 # Cost: 0.0912606 Accuracy: 100 / 100  
487 # Cost: 0.0909272 Accuracy: 100 / 100  
488 # Cost: 0.090496 Accuracy: 100 / 100  
489 # Cost: 0.0920291 Accuracy: 100 / 100  
490 # Cost: 0.0917305 Accuracy: 100 / 100  
491 # Cost: 0.0915103 Accuracy: 100 / 100  
492 # Cost: 0.0899408 Accuracy: 100 / 100  
493 # Cost: 0.0914209 Accuracy: 99 / 100  
494 # Cost: 0.092102 Accuracy: 100 / 100  
495 # Cost: 0.093399 Accuracy: 100 / 100  
496 # Cost: 0.0926048 Accuracy: 100 / 100  
497 # Cost: 0.0917633 Accuracy: 100 / 100  
498 # Cost: 0.0921034 Accuracy: 100 / 100  
499 # Cost: 0.0928722 Accuracy: 100 / 100  
500 # Cost: 0.0942567 Accuracy: 100 / 100  
501 # Cost: 0.091184 Accuracy: 100 / 100  
502 # Cost: 0.0913637 Accuracy: 100 / 100  
503 # Cost: 0.091654 Accuracy: 100 / 100  
504 # Cost: 0.0915468 Accuracy: 100 / 100  
505 # Cost: 0.0896783 Accuracy: 100 / 100  
506 # Cost: 0.0893992 Accuracy: 100 / 100  
507 # Cost: 0.0916488 Accuracy: 100 / 100  
508 # Cost: 0.0921669 Accuracy: 100 / 100  
509 # Cost: 0.0916404 Accuracy: 100 / 100  
510 # Cost: 0.0908221 Accuracy: 100 / 100  
511 # Cost: 0.0910142 Accuracy: 100 / 100  
512 # Cost: 0.0913771 Accuracy: 100 / 100  
513 # Cost: 0.09243 Accuracy: 100 / 100  
514 # Cost: 0.0912463 Accuracy: 100 / 100  
515 # Cost: 0.0911145 Accuracy: 100 / 100  
516 # Cost: 0.0911484 Accuracy: 100 / 100  
517 # Cost: 0.0907027 Accuracy: 100 / 100  
518 # Cost: 0.0911056 Accuracy: 100 / 100  
519 # Cost: 0.0914963 Accuracy: 100 / 100  
520 # Cost: 0.0926192 Accuracy: 100 / 100  
521 # Cost: 0.0907697 Accuracy: 99 / 100  
522 # Cost: 0.09215 Accuracy: 100 / 100  
523 # Cost: 0.0927574 Accuracy: 100 / 100  
524 # Cost: 0.0908841 Accuracy: 100 / 100  
525 # Cost: 0.0922874 Accuracy: 100 / 100  
526 # Cost: 0.0930886 Accuracy: 100 / 100  
527 # Cost: 0.092232 Accuracy: 100 / 100  
528 # Cost: 0.0916903 Accuracy: 100 / 100  
529 # Cost: 0.0924921 Accuracy: 100 / 100  
530 # Cost: 0.0908627 Accuracy: 100 / 100  
531 # Cost: 0.0904231 Accuracy: 100 / 100  
532 # Cost: 0.0914372 Accuracy: 100 / 100  
533 # Cost: 0.0906258 Accuracy: 100 / 100  
534 # Cost: 0.0898558 Accuracy: 100 / 100  
535 # Cost: 0.0919011 Accuracy: 100 / 100  
536 # Cost: 0.0915588 Accuracy: 100 / 100  
537 # Cost: 0.0921337 Accuracy: 100 / 100  
538 # Cost: 0.0908392 Accuracy: 100 / 100  
539 # Cost: 0.0925348 Accuracy: 100 / 100  
540 # Cost: 0.0915594 Accuracy: 100 / 100  
541 # Cost: 0.0906666 Accuracy: 100 / 100  
542 # Cost: 0.0916788 Accuracy: 100 / 100  
543 # Cost: 0.0915069 Accuracy: 100 / 100  
544 # Cost: 0.091643 Accuracy: 100 / 100  
545 # Cost: 0.0900198 Accuracy: 100 / 100  
546 # Cost: 0.0920953 Accuracy: 100 / 100  
547 # Cost: 0.0918621 Accuracy: 100 / 100  
548 # Cost: 0.0917129 Accuracy: 100 / 100  
549 # Cost: 0.0915582 Accuracy: 100 / 100  
550 # Cost: 0.092318 Accuracy: 100 / 100  
551 # Cost: 0.0909816 Accuracy: 100 / 100  
552 # Cost: 0.0922378 Accuracy: 100 / 100  
553 # Cost: 0.0908495 Accuracy: 100 / 100  
554 # Cost: 0.0914777 Accuracy: 100 / 100  
555 # Cost: 0.0907575 Accuracy: 100 / 100  
556 # Cost: 0.0918842 Accuracy: 100 / 100  
557 # Cost: 0.0913665 Accuracy: 100 / 100  
558 # Cost: 0.0903595 Accuracy: 100 / 100  
559 # Cost: 0.0912446 Accuracy: 100 / 100  
560 # Cost: 0.0910663 Accuracy: 100 / 100  
561 # Cost: 0.0909414 Accuracy: 100 / 100  
562 # Cost: 0.0916522 Accuracy: 100 / 100  
563 # Cost: 0.0907814 Accuracy: 100 / 100  
564 # Cost: 0.090983 Accuracy: 100 / 100  
565 # Cost: 0.0927666 Accuracy: 100 / 100  
566 # Cost: 0.090397 Accuracy: 100 / 100  
567 # Cost: 0.0917119 Accuracy: 100 / 100  
568 # Cost: 0.0898354 Accuracy: 100 / 100  
569 # Cost: 0.0918394 Accuracy: 100 / 100  
570 # Cost: 0.089002 Accuracy: 100 / 100  
571 # Cost: 0.0916902 Accuracy: 100 / 100  
572 # Cost: 0.092083 Accuracy: 100 / 100  
573 # Cost: 0.0915101 Accuracy: 100 / 100  
574 # Cost: 0.0931548 Accuracy: 100 / 100  
575 # Cost: 0.0893782 Accuracy: 100 / 100  
576 # Cost: 0.0913528 Accuracy: 100 / 100  
577 # Cost: 0.0921975 Accuracy: 100 / 100  
578 # Cost: 0.0894949 Accuracy: 100 / 100  
579 # Cost: 0.0936981 Accuracy: 100 / 100  
580 # Cost: 0.0919966 Accuracy: 100 / 100  
581 # Cost: 0.0892555 Accuracy: 100 / 100  
582 # Cost: 0.089168 Accuracy: 100 / 100  
583 # Cost: 0.0905726 Accuracy: 100 / 100  
584 # Cost: 0.0908773 Accuracy: 100 / 100  
585 # Cost: 0.0918277 Accuracy: 100 / 100  
586 # Cost: 0.0898744 Accuracy: 100 / 100  
587 # Cost: 0.0905348 Accuracy: 100 / 100  
588 # Cost: 0.0891129 Accuracy: 100 / 100  
589 # Cost: 0.0898332 Accuracy: 100 / 100  
590 # Cost: 0.0912857 Accuracy: 100 / 100  
591 # Cost: 0.0911505 Accuracy: 100 / 100  
592 # Cost: 0.0908127 Accuracy: 100 / 100  
593 # Cost: 0.0916399 Accuracy: 100 / 100  
594 # Cost: 0.0899792 Accuracy: 100 / 100  
595 # Cost: 0.0914398 Accuracy: 100 / 100  
596 # Cost: 0.089838 Accuracy: 100 / 100  
597 # Cost: 0.0895193 Accuracy: 100 / 100  
598 # Cost: 0.0909437 Accuracy: 100 / 100  
599 # Cost: 0.0926661 Accuracy: 100 / 100  
