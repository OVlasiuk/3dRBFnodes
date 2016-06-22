function points = irrational_nodes(corner1, corner2, num)

diagonal = corner2-corner1;
points = zeros(num,3);
for i=1:num
    points(i,:) = [diagonal(1)*i/num, diagonal(2)*frac_part(sqrt(2)*i), diagonal(3)*frac_part(sqrt(3) *i)];
    points(i,:) = points(i,:) + corner1;
end

% plot3(points(1,:), points(2,:), points(3,:),  '.k');