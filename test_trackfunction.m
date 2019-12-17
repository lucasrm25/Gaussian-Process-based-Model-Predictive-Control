[trackdata, x0, th0, w] = RaceTrack.loadTrack_02();

track = RaceTrack(trackdata, x0, th0, w);

track.plotTrack()
hold on;

N = max(track.dist)

[pos_c, psi_c, R_c] = track.getTrackInfo(1:1:N)
% plot(pos_c(1,:),pos_c(2,:),'-o')

figure
plot(rad2deg(psi_c))

scatter(track.track_c(1,:),track.track_c(2,:))



[pos_c, psi_c, R_c] = track.getTrackInfo(122)
% scatter(track.track_c(1,:),track.track_c(2,:))
scatter(pos_c(1,:),pos_c(2,:))


vehicle_pos = [0;-35];
vehicle_dist = track.getTrackDistance(vehicle_pos)

sc1 = scatter(vehicle_pos(1), vehicle_pos(2), 50, 'filled', 'DisplayName','vehicle pos')
[pos_c, psi_c, R_c] = track.getTrackInfo(vehicle_dist);
sc2 = scatter(pos_c(1),pos_c(2), 50, 'filled', 'DisplayName','closest track point')


figure
plot(track.dist)
ylim([122 123])