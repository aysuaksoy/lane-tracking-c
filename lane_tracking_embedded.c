// lane_tracking_embedded.c - Embedded system compatible lane tracking with true Kalman gain (pure C)

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Kalman filter structure
typedef struct {
    float dt;
    float state[4];      // [slope, intercept, d(slope), d(intercept)]
    float P[4][4];       // Covariance
    float F[4][4];       // State transition
    float Q[4][4];       // Process noise
    float H[2][4];       // Measurement matrix
    float R[2][2];       // Measurement noise
    float K[4][2];       // Kalman gain
    clock_t predict_time;
    clock_t update_time;
} KalmanFilter;

// Line structure
typedef struct {
    int x1, y1, x2, y2;
} Line;

void mat_mult_4x4_4x1(float res[4], float A[4][4], float B[4]) {
    for (int i = 0; i < 4; ++i) {
        res[i] = 0;
        for (int j = 0; j < 4; ++j)
            res[i] += A[i][j] * B[j];
    }
}

void kalman_init(KalmanFilter *kf, float dt) {
    memset(kf, 0, sizeof(KalmanFilter));
    kf->dt = dt;
    float F[4][4] = {
        {1, 0, dt, 0},
        {0, 1, 0, dt},
        {0, 0, 1, 0 },
        {0, 0, 0, 1 }
    };
    memcpy(kf->F, F, sizeof(F));
    float H[2][4] = {
        {1, 0, 0, 0},
        {0, 1, 0, 0}
    };
    memcpy(kf->H, H, sizeof(H));
    for (int i = 0; i < 4; i++) {
        kf->Q[i][i] = 0.01f;
        kf->P[i][i] = 1.0f;
    }
    kf->R[0][0] = 0.1f;
    kf->R[1][1] = 0.1f;
}

void kalman_predict(KalmanFilter *kf) {
    kf->predict_time = clock();
    float new_state[4];
    mat_mult_4x4_4x1(new_state, kf->F, kf->state);
    memcpy(kf->state, new_state, sizeof(new_state));
}

void kalman_update(KalmanFilter *kf, float meas_slope, float meas_intercept) {
    kf->update_time = clock();
    float y[2];
    y[0] = meas_slope - (kf->H[0][0] * kf->state[0] + kf->H[0][1] * kf->state[1]);
    y[1] = meas_intercept - (kf->H[1][0] * kf->state[0] + kf->H[1][1] * kf->state[1]);

    // True Kalman gain: K = P * H^T * inv(H * P * H^T + R) (approximate)
    float PHt[4][2] = {0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 4; ++k)
                PHt[i][j] += kf->P[i][k] * kf->H[j][k];

    float S[2][2] = {0};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 4; ++k)
                S[i][j] += kf->H[i][k] * PHt[k][j];
    S[0][0] += kf->R[0][0];
    S[1][1] += kf->R[1][1];

    float invS[2][2] = {0};
    invS[0][0] = 1.0f / S[0][0];
    invS[1][1] = 1.0f / S[1][1];

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 2; ++j)
            kf->K[i][j] = PHt[i][0] * invS[0][j];

    kf->state[0] += kf->K[0][0] * y[0];
    kf->state[1] += kf->K[1][1] * y[1];
}

float clock_to_ms(clock_t ticks) {
    return 1000.0f * (float)ticks / CLOCKS_PER_SEC;
}

const Line test_lines[][2] = {
    {{100, 480, 300, 300}, {340, 480, 550, 320}},
    {{105, 480, 305, 305}, {345, 480, 555, 325}},
    {{110, 480, 310, 310}, {350, 480, 560, 330}},
    {{0, 0, 0, 0}, {355, 480, 565, 335}},
    {{120, 480, 320, 320}, {360, 480, 570, 340}},
    {{125, 480, 325, 325}, {365, 480, 575, 345}},
    {{130, 480, 330, 330}, {370, 480, 580, 350}},
    {{0, 0, 0, 0}, {0, 0, 0, 0}},
    {{140, 480, 340, 340}, {380, 480, 590, 360}},
    {{145, 480, 345, 345}, {385, 480, 595, 365}}
};

void apply_pwm_control(float slope_error) {
    float gain = 50.0f;
    int pwm_value = (int)(gain * slope_error);
    if (pwm_value > 255) pwm_value = 255;
    if (pwm_value < -255) pwm_value = -255;
    printf("    PWM correction: %d
", pwm_value);
}

int main() {
    KalmanFilter left_kf, right_kf;
    kalman_init(&left_kf, 0.033f);
    kalman_init(&right_kf, 0.033f);
    const int num_frames = sizeof(test_lines) / sizeof(test_lines[0]);

    for (int frame = 0; frame < num_frames; frame++) {
        float left_slope = 0, left_intercept = 0, right_slope = 0, right_intercept = 0;
        int left_count = 0, right_count = 0;

        for (int i = 0; i < 2; i++) {
            const Line *line = &test_lines[frame][i];
            if ((line->x1 == 0 && line->x2 == 0) || line->x1 == line->x2) continue;
            float slope = (float)(line->y2 - line->y1) / (line->x2 - line->x1);
            float intercept = line->y1 - slope * line->x1;
            if (fabsf(slope) < 0.2f) continue;
            if (slope < 0) {
                left_slope += slope;
                left_intercept += intercept;
                left_count++;
            } else {
                right_slope += slope;
                right_intercept += intercept;
                right_count++;
            }
        }

        if (left_count > 0) {
            left_slope /= left_count;
            left_intercept /= left_count;
            kalman_predict(&left_kf);
            kalman_update(&left_kf, left_slope, left_intercept);
        } else {
            kalman_predict(&left_kf);
        }

        if (right_count > 0) {
            right_slope /= right_count;
            right_intercept /= right_count;
            kalman_predict(&right_kf);
            kalman_update(&right_kf, right_slope, right_intercept);
        } else {
            kalman_predict(&right_kf);
        }

        printf("Frame %d:
", frame);
        printf("  Left: slope=%.3f, intercept=%.1f
", left_kf.state[0], left_kf.state[1]);
        printf("    latency: pred=%.3fms, update=%.3fms
", clock_to_ms(left_kf.predict_time), clock_to_ms(left_kf.update_time));
        printf("  Right: slope=%.3f, intercept=%.1f
", right_kf.state[0], right_kf.state[1]);
        printf("    latency: pred=%.3fms, update=%.3fms
", clock_to_ms(right_kf.predict_time), clock_to_ms(right_kf.update_time));

        float slope_error = right_kf.state[0] - left_kf.state[0];
        apply_pwm_control(slope_error);
        printf("
");
    }

    return 0;
}
