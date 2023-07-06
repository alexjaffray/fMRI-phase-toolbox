function plotSphericalHarmonics(shmat,time)

% Plot Each Coefficient of Varying Order
figure()
subplot(4,1,1);
plot(time,shmat(:,1));
legend('1');
subplot(4,1,2);
plot(time,shmat(:,2:4));
legend('x','y','z');
subplot(4,1,3);
plot(time,shmat(:,5:9));
legend('xy','yz','xz','x^2-y^2','2z^2-x^2-y^2');
subplot(4,1,4);
plot(time,shmat(:,10:16));
legend('xyz','zx^2-zy^2','3yx^2-y^3','(5z^2-(x^2+y^2+z^2))y','(5z^2-(x^2+y^2+z^2))x','5z^3-3(x^2+y^2+z^2)z','x^3-3xy^2');

end

