// HCSR04wrapper.cpp  (Non-scanning, NewPing-free, pulseIn-based)
#include <Arduino.h>

static int g_trig_pin[3] = { -1, -1, -1 };
static int g_echo_pin[3] = { -1, -1, -1 };
static unsigned long g_timeout_us[3] = { 12000UL, 12000UL, 12000UL };

extern "C" void HCSR04Sonar_Init(int trigger_pin, int echo_pin, int Sonar)
{
    int idx = Sonar - 1;
    if (idx < 0 || idx >= 3) return;

    g_trig_pin[idx] = trigger_pin;
    g_echo_pin[idx] = echo_pin;

    pinMode(trigger_pin, OUTPUT);
    digitalWrite(trigger_pin, LOW);
    pinMode(echo_pin, INPUT); 
}

static inline void trigger_pulse(uint8_t trig) {
    digitalWrite(trig, LOW);
    delayMicroseconds(2);
    digitalWrite(trig, HIGH);
    delayMicroseconds(10);
    digitalWrite(trig, LOW);
}

extern "C" int HCSR04Sonar_Read(int Sonar)
{
    int idx = Sonar - 1;
    if (idx < 0 || idx >= 3) return 0;
    int trig = g_trig_pin[idx], echo = g_echo_pin[idx];
    if (trig < 0 || echo < 0) return 0;

    trigger_pulse((uint8_t)trig);

    unsigned long dur = pulseIn((uint8_t)echo, HIGH, g_timeout_us[idx]);
    if (dur == 0UL) return 0;

    // int cm = (int)(dur / 58UL);
    int cm = (int)(dur);
    if (cm < 0) cm = 0;
    if (cm > 500) cm = 500;
    return cm;
}
