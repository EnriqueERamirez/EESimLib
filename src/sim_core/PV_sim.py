import numpy as np
import math
from enum import Enum
from datetime import datetime, timedelta

class TecnologiaPanel(Enum):
    """Tipos de tecnologías de paneles solares con sus coeficientes de temperatura"""
    SILICIO_CRISTALINO = -0.0035
    CIS = -0.0030
    CDTE = -0.0025

class ParametrosElectricos:
    """Clase para almacenar parámetros eléctricos básicos"""
    def __init__(self, tension, intensidad, potencia):
        self.tension = tension
        self.intensidad = intensidad
        self.potencia = potencia

    def __str__(self):
        return (f"Tensión: {self.tension:.2f}V\n"
                f"Intensidad: {self.intensidad:.2f}A\n"
                f"Potencia: {self.potencia:.2f}W")

class ParametrosElectricosIrradiancia(ParametrosElectricos):
    """Extiende ParametrosElectricos para incluir irradiancia"""
    def __init__(self, tension, intensidad, potencia, irradiancia):
        super().__init__(tension, intensidad, potencia)
        self.irradiancia = irradiancia

    def __str__(self):
        return (f"Irradiancia: {self.irradiancia} W/m²\n" + super().__str__())

class PotenciaPanelFV:
    """Clase para almacenar resultados de potencia y eficiencia del panel"""
    def __init__(self, potencia_dc, eficiencia, irradiancia, temperatura):
        self.potencia_dc = potencia_dc
        self.eficiencia = eficiencia
        self.irradiancia = irradiancia
        self.temperatura = temperatura

    def __str__(self):
        return (f"Temperatura: {self.temperatura}°C\n"
                f"Irradiancia: {self.irradiancia} W/m²\n"
                f"Eficiencia: {self.eficiencia*100:.2f}%\n"
                f"Potencia DC: {self.potencia_dc:.2f}W")

class PanelSolar:
    """Clase principal para representar un panel solar"""
    def __init__(self,
                 area,                  # Área del módulo en m²
                 eficiencia_stc,        # Eficiencia en condiciones STC
                 V_oc_stc,              # Tensión circuito abierto STC
                 I_sc_stc,              # Intensidad cortocircuito STC
                 P_pmp_stc,             # Potencia máxima STC
                 FF,                    # Factor de forma
                 tonc=47,               # Temperatura de operación nominal
                 tecnologia=TecnologiaPanel.SILICIO_CRISTALINO):
        self.area = area
        self.eficiencia_stc = eficiencia_stc
        self.V_oc_stc = V_oc_stc
        self.I_sc_stc = I_sc_stc
        self.P_pmp_stc = P_pmp_stc
        self.FF = FF
        self.tonc = tonc
        self.tecnologia = tecnologia

class ResultadoSimulacion:
    """Clase para almacenar resultados de la simulación temporal"""
    def __init__(self, tiempo, temp_ambiente, irradiancia, temp_celula,
                 tension, intensidad, potencia, eficiencia):
        self.tiempo = tiempo
        self.temp_ambiente = temp_ambiente
        self.irradiancia = irradiancia
        self.temp_celula = temp_celula
        self.tension = tension
        self.intensidad = intensidad
        self.potencia = potencia
        self.eficiencia = eficiencia

    def resumen_diario(self):
        """Genera un resumen de los resultados diarios"""
        idx_max_power = np.argmax(self.potencia)
        energia_total = np.trapz(self.potencia, dx=1/6) / 1000  # kWh (para intervalos de 10 min)

        return {
            'potencia_maxima': self.potencia[idx_max_power],
            'hora_maxima_potencia': self.tiempo[idx_max_power],
            'temp_max_celula': np.max(self.temp_celula),
            'eficiencia_promedio': np.mean(self.eficiencia) * 100,
            'energia_total': energia_total
        }

def calcular_temperatura_celula(temp_ambiente, irradiancia, tonc=47):
    """
    Calcula la temperatura de la célula fotovoltaica.

    Args:
        temp_ambiente (float): Temperatura ambiente en °C
        irradiancia (float): Irradiancia solar en W/m²
        tonc (float): Temperatura de operación nominal de la célula

    Returns:
        float: Temperatura de la célula en °C
    """
    return temp_ambiente + irradiancia * (tonc - 20) / 800

def calcular_parametros_temperatura(V_oc, I_sc, P_pmp, T_m, T_stc=25,
                                  alpha=0.0004, beta=-0.0025, gamma=-0.0037):
    """
    Calcula parámetros eléctricos considerando la temperatura.
    """
    delta_T = T_m - T_stc
    V_oc_tm = V_oc * (1 + beta * delta_T)
    I_sc_tm = I_sc * (1 + alpha * delta_T)
    P_pmp_tm = P_pmp * (1 + gamma * delta_T)

    return ParametrosElectricos(V_oc_tm, I_sc_tm, P_pmp_tm)

def calcular_parametros_irradiancia(V_oc_stc, I_sc_stc, FF, G_m, T=298.15,
                                  G_stc=1000, n=1.2):
    """
    Calcula parámetros eléctricos considerando la irradiancia.
    Maneja casos de baja irradiancia o nocturnos.
    """
    k_b = 1.381e-23  # Constante de Boltzmann
    q = 1.602e-19    # Carga del electrón

    # Manejar casos de irradiancia muy baja o nula
    if G_m <= 1:  # Consideramos irradiancia mínima de 1 W/m²
        return ParametrosElectricosIrradiancia(0, 0, 0, G_m)

    # Cálculos normales para irradiancia suficiente
    V_oc_gm = V_oc_stc + (n * k_b * T / q) * math.log(G_m / G_stc)
    I_sc_gm = I_sc_stc * (G_m / G_stc)
    P_pmp_gm = FF * V_oc_gm * I_sc_gm

    return ParametrosElectricosIrradiancia(V_oc_gm, I_sc_gm, P_pmp_gm, G_m)

def calcular_potencia_panel(eficiencia_stc, area, G_m, T_m,
                          tecnologia=TecnologiaPanel.SILICIO_CRISTALINO,
                          T_stc=25):
    """
    Calcula la potencia de salida considerando temperatura e irradiancia.
    """
    k = tecnologia.value
    eficiencia_temp = eficiencia_stc * (1 + k * (T_m - T_stc))
    potencia_dc = eficiencia_temp * area * G_m

    return PotenciaPanelFV(potencia_dc, eficiencia_temp, G_m, T_m)

def simular_panel_solar(panel: PanelSolar,
                       temp_ambiente: np.ndarray,
                       irradiancia: np.ndarray,
                       tiempo: np.ndarray) -> ResultadoSimulacion:
    """
    Simula el comportamiento temporal del panel solar.
    """
    n_points = len(tiempo)
    temp_celula = np.zeros(n_points)
    tension = np.zeros(n_points)
    intensidad = np.zeros(n_points)
    potencia = np.zeros(n_points)
    eficiencia = np.zeros(n_points)

    for i in range(n_points):
        # Temperatura de la célula
        temp_celula[i] = calcular_temperatura_celula(
            temp_ambiente[i],
            irradiancia[i],
            panel.tonc
        )

        # Efectos de la temperatura
        params_temp = calcular_parametros_temperatura(
            panel.V_oc_stc,
            panel.I_sc_stc,
            panel.P_pmp_stc,
            temp_celula[i]
        )

        # Efectos de la irradiancia
        params_irrad = calcular_parametros_irradiancia(
            panel.V_oc_stc,
            panel.I_sc_stc,
            panel.FF,
            irradiancia[i]
        )

        # Potencia y eficiencia final
        resultado_potencia = calcular_potencia_panel(
            panel.eficiencia_stc,
            panel.area,
            irradiancia[i],
            temp_celula[i],
            panel.tecnologia
        )

        # Almacenar resultados
        tension[i] = params_irrad.tension
        intensidad[i] = params_irrad.intensidad
        potencia[i] = resultado_potencia.potencia_dc
        eficiencia[i] = resultado_potencia.eficiencia

    return ResultadoSimulacion(
        tiempo, temp_ambiente, irradiancia, temp_celula,
        tension, intensidad, potencia, eficiencia
    )

def generar_perfil_dia(fecha_base=None, intervalo_minutos=10):
    """
    Genera perfiles diarios de temperatura e irradiancia.
    """
    if fecha_base is None:
        fecha_base = datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)

    # Generar timestamps
    n_points = int(24 * 60 / intervalo_minutos)
    horas = np.linspace(0, 24, n_points)
    tiempo = np.array([fecha_base + timedelta(hours=h) for h in horas])

    # Generar temperatura (patrón sinusoidal)
    temp_ambiente = 20 + 10 * np.sin(np.pi * (horas - 6) / 12)

    # Generar irradiancia (patrón gaussiano)
    irradiancia = 1000 * np.exp(-((horas - 12)**2) / 20)
    irradiancia = np.maximum(irradiancia, 0)  # No permitir valores negativos

    return tiempo, temp_ambiente, irradiancia
