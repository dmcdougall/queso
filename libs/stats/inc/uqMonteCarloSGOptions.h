/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_MOC_SG_OPTIONS_H__
#define __UQ_MOC_SG_OPTIONS_H__

#include <uqEnvironment.h>
#include <uqSequenceStatisticalOptions.h>

class uqMonteCarloSGOptionsClass
{
public:
  uqMonteCarloSGOptionsClass(const uqBaseEnvironmentClass& env, const char* prefix);
 ~uqMonteCarloSGOptionsClass();

  void scanOptionsValues();
  void print            (std::ostream& os) const;

  std::string                        m_prefix;

private:
  void   defineMyOptions  (po::options_description& optionsDesc) const;
  void   getMyOptionValues(po::options_description& optionsDesc);

  const uqBaseEnvironmentClass& m_env;
  po::options_description*      m_optionsDesc;

  std::string                   m_option_help;
};

std::ostream& operator<<(std::ostream& os, const uqMonteCarloSGOptionsClass& obj);
#endif // __UQ_MOC_SG_OPTIONS_H__
