% 1. 使用推荐的 datetime 函数获取当前时间，并设置格式
currentTime = datetime('now', 'Format', 'MM/dd HH:mm'); % 格式化为 "月/日 时:分"
% 转换为字符向量用于命令
timestamp = char(currentTime);

% 2. 构建完整的Git命令并执行
commit_msg = ['Update MATLAB scripts ', timestamp];
command = ['git add *.m && git commit -m "', commit_msg, '" && git push'];
system(command);

% 3. 可选：显示执行结果
fprintf('已执行命令: %s\n', command);