function xydistribution(app,x,y)
c = ksdensity([x,y], [x,y]);
cmp = [208	209	230
    166	189	219
    116	169	207
    54	144	192
    5	112	176
    4	90	141
    2	56	88
    ]/255;
colormap(app.hs,cmp)
Gs = round(interp1(linspace(min(c),max(c),size(cmp,1)),1:size(cmp,1),c));
colors = cmp(Gs,:); % make RGB matrix from scaled.

scatter(app.hs,x, y, 30, colors,'filled');
histogram(app.hx,x);
histogram(app.hy,y,'Orientation','horizontal');

box(app.hs,'on')
box(app.hx,'off')
box(app.hy,'off')
xticklabels(app.hx,{})
xticklabels(app.hy,{})
yticklabels(app.hx,{})
yticklabels(app.hy,{})
app.hy.YColor='k';
app.hx.XColor='k';
app.hy.XColor='none';
app.hx.YColor='none';